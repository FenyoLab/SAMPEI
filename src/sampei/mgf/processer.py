import time

from pyteomics import mgf

from src.sampei.mgf.MGF import MGF


def process(mgf_path, num_intensities=20):
    """Load the mgf file into an array

    :param mgf_path: The path of the mgf file to be loaded
    :returns: A list of mgf entries

    """
    MGF_FILE = mgf.read(mgf_path)
    entries = []
    count = 0
    for f in MGF_FILE:
        tmp = MGF(f["params"], f["m/z array"], f["intensity array"], num_intensities)
        entries.append(tmp)
        count += 1
    # Sort the mgf entries using the formula used for lookups to allow for binary searching
    return sorted(entries, key=lambda mgf: mgf.pepmass_charge)


def main(path):
    process(path)


if __name__ == "__main__":
    main("./data/mgf/171202_Ecoli_ctrl2_1ug.mgf")
