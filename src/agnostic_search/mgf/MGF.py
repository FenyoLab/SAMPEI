import math
import numpy as np

from collections import namedtuple
from typing import Dict

OVERLAP = namedtuple(
    "Overlap",
    [
        "t_scan",
        "q_scan",
        "count",
        "target_matched",
        "query_matched",
        "target_log_sum",
        "query_log_sum",
        "difference",
    ],
)


class MGF:
    __slots__ = [
        "title",
        "scan",
        "rts",
        "pepmass",
        "pepmass_intensity",
        "pepmass_charge",
        "charge",
        "intensities",
        "m_zs",
        "mz_sum",
        "intensities_sum",
        "top_idx",
    ]

    def __init__(
        self,
        params: Dict[str, str],
        mass_per_charges: np.ndarray,
        intensities: np.ndarray,
        num_intensities: int = 20,
    ):
        self.title = params["title"]
        # FIXME: Scan is not always desegnated using .
        self.scan = int(self.title.split(".")[1])
        self.rts = params["rtinseconds"]
        self.pepmass = params["pepmass"][0]
        self.pepmass_intensity = (
            params["pepmass"][1] if len(params["pepmass"]) > 1 else 0
        )
        # self.charge = params["charge"][0].denominator
        self.charge = params["charge"][0].real
        self.m_zs = mass_per_charges
        self.intensities = intensities

        self.pepmass_charge = self.pepmass * self.charge

        self.top_idx = self.get_top_intensity_indices(num_intensities)[::-1]

        self.mz_sum = self.m_zs[self.top_idx].sum()
        self.intensities_sum = self.intensities[self.top_idx].sum()

    def get_top_intensity_indices(self, n: int = 20) -> np.ndarray:
        top_idx = self.intensities.argsort()[-n:]
        return top_idx

    def get_window(self, window: int = 200):
        mz_range = window / self.charge
        return (
            (self.pepmass - mz_range) * self.charge,
            (self.pepmass + mz_range) * self.charge,
        )

    # TODO: Brute force approach needs to be optimized
    def count_matches(self, query, E, measure):
        if measure == "ppm":
            ppm = 1e6
        difference = self.pepmass_charge - query.pepmass_charge
        count = 0
        target_total = 0
        query_total = 0
        # FIXME: This algorithm should improve the performance
        #      : ALTERNATIVELY a binary sort would work
        #     do {
        #     status = compare(shortList[i], longList[j]);
        #     if(status == EQUAL) {
        #         # Found match!
        #         i++;
        #         j++;
        #     } else if(status == EARLIER) {
        #         # No match, discard first entry in short list
        #         i++;
        #     } else {
        #         # No match, discard first entry in long list
        #         j++;
        #     }
        # } while( (i < shortListEntries) && (j < longListEntries) );
        if self.scan == 45992:
            print(self.intensities)
            print(self.top_idx)
            print(self.intensities[self.top_idx])
            # FIXME: THIS IS WHY I GET DIFFERENT RESULTS
            print(difference)
        for query_idx in query.top_idx:
            for target_idx in self.top_idx:

                if (
                    abs(self.m_zs[target_idx] - query.m_zs[query_idx])
                    < (E) / ppm * query.m_zs[query_idx]
                    or abs(self.m_zs[target_idx] - query.m_zs[query_idx] + difference)
                    < (E) / ppm * query.m_zs[query_idx]
                ):
                    count += 1
                    target_total += self.intensities[target_idx]
                    query_total += query.intensities[query_idx]
                    break
        return OVERLAP(
            self.scan,
            query.scan,
            count,
            target_total / self.intensities_sum,
            query_total / query.intensities_sum,
            np.log(self.intensities_sum) / np.log(10),
            np.log(query.intensities_sum) / np.log(10),
            difference,
        )

    def __str__(self):
        return (
            "Title: {0.title}\n"
            "\tScan: {0.scan}\n"
            "\tCharge: {0.charge}\n"
            "\tRtins: {0.rts}\n"
            "\tPepmass: {0.pepmass}\n"
        ).format(self)
