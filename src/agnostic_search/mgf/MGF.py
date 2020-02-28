import numpy as np

from typing import Dict


class MGF:
    __slots__ = [
        "title",
        "scan",
        "rts",
        "pepmass",
        "pepmass_intensity",
        "charge",
        "intensities",
        "m_zs",
    ]

    def __init__(
        self,
        params: Dict[str, str],
        mass_per_charges: np.ndarray,
        intensities: np.ndarray,
    ):
        self.title = params["title"]
        # FIXME: Scan is not always desegnated using .
        self.scan = self.title.split(".")[1]
        self.rts = params["rtinseconds"]
        self.pepmass = params["pepmass"][0]
        self.pepmass_intensity = (
            params["pepmass"][1] if len(params["pepmass"]) > 1 else 0
        )
        # self.charge = params["charge"][0].denominator
        self.charge = params["charge"][0].real
        self.m_zs = mass_per_charges
        self.intensities = intensities

    def get_top_intensities(self, n: int = 20) -> np.ndarray:
        top_idx = self.intensities.argsort()[-n:]
        return self.intensities[top_idx]

    def get_window(self, window: int = 200):
        mz_range = window / self.charge
        return (
            (self.pepmass - mz_range) * self.charge,
            (self.pepmass + mz_range) * self.charge,
        )

    def __str__(self):
        return (
            "Title: {0.title}\n"
            "\tScan: {0.scan}\n"
            "\tCharge: {0.charge}\n"
            "\tRtins: {0.rts}\n"
            "\tPepmass: {0.pepmass}\n"
        ).format(self)
