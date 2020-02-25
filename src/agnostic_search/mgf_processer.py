from pyteomics import mgf


def main():
    mgf_path = "/home/mark/Data1/agnostic_search/mgf/171202_Ecoli_ctrl2_1ug.mgf"
    MGF_FILE = mgf.read(mgf_path)
    for f in MGF_FILE:
        print(type(f))
        break


if __name__ == "__main__":
    main()
