from pyteomics import mgf


def main(mgf_path):
    MGF_FILE = mgf.read(mgf_path)
    for f in MGF_FILE:
        print(f["params"]["title"])


if __name__ == "__main__":
    main("./data/mgf/171202_Ecoli_ctrl2_1ug.mgf")
