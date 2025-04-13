import math
import numpy as np
import re
import time

from numba import jit

from collections import namedtuple
from typing import Dict


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print("%r  %2.2f s" % (method.__name__, (te - ts)))
        return result

    return timed


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
        "time",
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

        _scan_match = re.findall(r'Scan\W(\S*?),', self.title)
        if _scan_match:
            self.scan = int(_scan_match[0])
            _time_match = re.findall(r'Time=(\S*?),', self.title)
            self.time = _time_match[0] if _time_match else 'NA'
        else:
            self.scan = int(self.title.split(".")[1])
            self.time = 'NA'
        self.rts = 'NA'
        if "rtinseconds" in params:
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

    #Update to different range (f.e. -50 to +500 Da)
    def get_window(self, upper_window: int = 200, lower_window: int = -50):
        return (
            (self.pepmass + lower_window) * self.charge,
            (self.pepmass + upper_window) * self.charge,
        )

    @staticmethod
    @jit(nopython=True)
    def count_numba(t_mzs, q_mzs, t_int, q_int, tidxs, qidxs, difference, E, ppm):
        count = 0
        target_total = 0
        query_total = 0
        for query_idx in qidxs:
            for target_idx in tidxs:

                if (
                    abs(t_mzs[target_idx] - q_mzs[query_idx])
                    < (E) / ppm * q_mzs[query_idx]
                    or abs(t_mzs[target_idx] - (q_mzs[query_idx] + difference))
                    < (E) / ppm * q_mzs[query_idx]
                ):
                    count += 1
                    target_total += t_int[target_idx]
                    query_total += q_int[query_idx]
                    break

        return (
            count,
            target_total,
            query_total,
        )

    def count_matches(self, query, E, measure):
        if measure == "ppm":
            ppm = 1e6
        difference = self.pepmass_charge - query.pepmass_charge
        count, target_total, query_total = self.count_numba(
            self.m_zs,
            query.m_zs,
            self.intensities,
            query.intensities,
            self.top_idx,
            query.top_idx,
            difference,
            E,
            ppm,
        )

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
            "\tPepmass Charge: {0.pepmass_charge}\n"
        ).format(self)
