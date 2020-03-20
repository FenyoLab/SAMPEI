import numpy as np
import re
import time

from pyteomics import mgf

from MGF import MGF
from processer import process


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print("%r  %2.2f s" % (method.__name__, (te - ts)))
        return result

    return timed


LINE_MAP = {
    "TITLE": lambda x: x,
    "RTINSECONDS": float,
    "PEPMASS": lambda x: tuple(float(y) for y in x.split()),
    "CHARGE": lambda x: int(x[:-1]),
}


@timeit
def parse_python(path):
    result = []
    scan = []
    with open(path) as f:
        block = False
        for line in f.readlines():
            line = line.strip()
            if line == "BEGIN IONS":
                block = True
                scan = {}
                mzs = []
                continue
            if block:
                if line == "END IONS":
                    scan["mzs"] = mzs
                    result.append(scan)
                    block = False
                    continue
                header = line.split("=")
                if header[0] in LINE_MAP:
                    scan[header[0]] = LINE_MAP[header[0]](header[1])
                else:
                    m_z, intensity = header[0].split()
                    mzs.append((float(m_z), float(intensity)))
    print(len(result))
    # print(result[0])


@timeit
def parse_numpy(path):
    regexp = re.compile(r"(BEGIN IONS.+?END IONS)+", re.MULTILINE | re.DOTALL)
    a = np.fromregex(path, regexp, [("match", "S2000")])
    print(len(a))


REGEX_MAP = {
    "TITLE": lambda x: re.findall(r"TITLE=(.+)", x),
    "RTINSECONDS": lambda x: re.findall(r"RTINSECONDS=(\d+\.\d+)", x),
    "PEPMASS": lambda x: re.findall(r"PEPMASS=(\d+\.\d+) (\d+\.\d+)", x),
    "CHARGE": lambda x: re.findall(r"CHARGE=(\d)", x),
    "mzi": lambda x: re.findall(
        re.compile(r"^(\d+?\.\d+) (\d+?\.\d+)$", re.MULTILINE), x
    ),
}


@timeit
def parse_regex(path):
    regexp = re.compile(r"(BEGIN IONS.+?END IONS)", re.MULTILINE | re.DOTALL)
    with open(path) as f:
        a = [
            {key: REGEX_MAP[key](m.group(0)) for key in REGEX_MAP}
            for m in regexp.finditer(f.read())
        ]
    print(len(a))
    # print(a[0])
    # print(re_title.find(a[0]))
    # print(re_mz_i.findall(a[0]))


@timeit
def parse_pyteomics(path):
    MGF_FILE = mgf.read(path)
    print(len(MGF_FILE))
    print(MGF_FILE[0])
    print(MGF_FILE[-1])


@timeit
def parse_mgf(path):
    MGF_FILE = mgf.read(path)
    entries = []
    count = 0
    for f in MGF_FILE:
        tmp = MGF(f["params"], f["m/z array"], f["intensity array"],)
        # TODO: Num top intensities should be a parameter
        entries.append(tmp)
        # if count == 5:
        #     break
        count += 1
    # sorted(entries, key=lambda mgf: mgf.pepmass * mgf.charge)
    return sorted(entries, key=lambda mgf: mgf.pepmass_charge)


def main(path):
    # parse_python(path)
    # parse_numpy(path)
    # parse_regex(path)
    # parse_pyteomics(path)
    print(parse_mgf(path)[0])
    print(process(path)[0])

    # regexp = re.compile(r"(BEGIN IONS.+?END IONS)", re.MULTILINE | re.DOTALL)
    # with open(path) as f:
    #     a = [m.group(1) for m in regexp.finditer(f.read())]
    # re_title = re.compile(r"^TITLE=(.+)$", re.MULTILINE)
    # re_mz_i = re.compile(r"^(\d+?\.\d+) (\d+?\.\d+)$", re.MULTILINE)
    # print(a[0])
    # print(re_title.find(a[0]))
    # print(re_mz_i.findall(a[0]))
    # regexp = re.compile(r"(BEGIN IONS.+?(\d+)END IONS)+", re.MULTILINE | re.DOTALL)
    # a = np.fromregex(path, regexp, [("match", "S2000"), ("ma", np.int32)])
    # print(a[0])


if __name__ == "__main__":
    PATH = "./data/mgf/171202_Ecoli_ctrl2_1ug.mgf"
    main(PATH)

