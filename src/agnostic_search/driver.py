import math
import numpy as np
import os
import pandas as pd
import pprint
import time

from src.agnostic_search.mgf.processer import process as mgf_processor
from src.agnostic_search.masses.masses import (
    MASSES,
    calc_peptide_mass,
    calc_peptide_fragments,
)


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print("%r  %2.2f s" % (method.__name__, (te - ts)))
        return result

    return timed


def find_query_range(mgf_queries, low, high):
    lower_bound = 0
    upper_bound = len(mgf_queries)
    while lower_bound < upper_bound:
        mid = lower_bound + (upper_bound - lower_bound) // 2
        if low <= mgf_queries[mid].pepmass * mgf_queries[mid].charge:
            upper_bound = mid
        else:
            lower_bound = mid + 1
    low_idx = lower_bound
    lower_bound = 0
    upper_bound = len(mgf_queries)

    while lower_bound < upper_bound:
        mid = lower_bound + (upper_bound - lower_bound) // 2
        if high < mgf_queries[mid].pepmass * mgf_queries[mid].charge:
            upper_bound = mid
        else:
            lower_bound = mid + 1
    # Add 1 so that the last index is inclusive when using slice
    high_idx = upper_bound

    return low_idx, high_idx


def calculate_largest_gap(pep, sequence_evidence_b, sequence_evidence_y):
    largest_gap = 0
    gap = 0
    for k in range(len(pep) - 1):
        if sequence_evidence_b[k] > 0 or sequence_evidence_y[len(pep) - 2 - k] > 0:
            gap = 0
        else:
            gap += 1
            if largest_gap < gap:
                largest_gap = gap
        return largest_gap


def select_fragments_scan(
    scan,
    frag_mz_min,
    mgf_pepmatch_ions,
    sequence_evidence_b,
    sequence_evidence_y,
    fragments,
    precursor_error,
    precursor_error_type,
    fragment_error,
    fragment_error_type,
    num_intense,
    mgf_table,
    current_target,
):
    matched_intensity = 0.0
    matched_intensity_not_parent = 0.0
    matched_intensity_parent = 0.0
    # ions = mgf_table[scan]
    ions = list(
        zip(
            current_target.m_zs[current_target.top_idx],
            current_target.intensities[current_target.top_idx],
        )
    )
    # sort ions for getting most intense
    if num_intense > 0:
        ions.sort(
            key=lambda tup: tup[1], reverse=True
        )  # sorts descending, by intensity
    intense_ions = ions[:num_intense]
    # sort ions for quicker compare
    ions.sort(key=lambda tup: tup[0])  # sorts ascending, by mz
    used_ions = np.zeros(len(ions))
    for frag_mz in fragments:
        error = fragment_error
        error_type = fragment_error_type
        if frag_mz[1].lstrip("-").startswith("[M+H"):
            error = precursor_error
            error_type = precursor_error_type
        count = 0
        for ion in ions:
            if error_type == "ppm":
                if abs(ion[0] - frag_mz[0]) < error * frag_mz[0] / 1e6:
                    mgf_pepmatch_ions.append([frag_mz[1], frag_mz[0], ion[0], ion[1]])
                    if frag_mz[1].lstrip("-").startswith("b"):
                        sequence_evidence_b[
                            int(
                                frag_mz[1]
                                .split("+")[0]
                                .split("o")[0]
                                .split("i")[0]
                                .split("*")[0]
                                .lstrip("-")
                                .lstrip("b")
                            )
                            - 1
                        ] += ion[1]
                    if frag_mz[1].lstrip("-").startswith("y"):
                        sequence_evidence_y[
                            int(
                                frag_mz[1]
                                .split("+")[0]
                                .split("o")[0]
                                .split("i")[0]
                                .split("*")[0]
                                .lstrip("-")
                                .lstrip("y")
                            )
                            - 1
                        ] += ion[1]
                    if used_ions[count] == 0 and frag_mz_min < ion[0]:
                        matched_intensity += ion[1]
                        if frag_mz[1].lstrip("-").startswith("["):
                            matched_intensity_parent += ion[1]
                        else:
                            matched_intensity_not_parent += ion[1]
                    used_ions[count] = 1
                elif frag_mz[0] < ion[0]:
                    break
            else:
                if abs(ion[0] - frag_mz[0]) < error:
                    mgf_pepmatch_ions.append([frag_mz[1], frag_mz[0], ion[0], ion[1]])
                    if frag_mz[1].lstrip("-").startswith("b"):
                        sequence_evidence_b[
                            int(
                                frag_mz[1]
                                .split("+")[0]
                                .split("o")[0]
                                .split("i")[0]
                                .split("*")[0]
                                .lstrip("-")
                                .lstrip("b")
                            )
                            - 1
                        ] += ion[1]
                    if frag_mz[1].lstrip("-").startswith("y"):
                        sequence_evidence_y[
                            int(
                                frag_mz[1]
                                .split("+")[0]
                                .split("o")[0]
                                .split("i")[0]
                                .split("*")[0]
                                .lstrip("-")
                                .lstrip("y")
                            )
                            - 1
                        ] += ion[1]
                    if used_ions[count] == 0 and frag_mz_min < ion[0]:
                        matched_intensity += ion[1]
                        if frag_mz[1].lstrip("-").startswith("["):
                            matched_intensity_parent += ion[1]
                        else:
                            matched_intensity_not_parent += ion[1]
                    used_ions[count] = 1
                elif frag_mz[0] < ion[0]:
                    break
            count += 1
    if 1 == 1:  # fragment_error_type=='ppm':
        mgf_pepmatch_ions_isotope1 = []
        for ion_type, mz_calc, mz_exp, intensity in mgf_pepmatch_ions:
            error = fragment_error
            error_type = fragment_error_type
            if ion_type.lstrip("-").startswith("[M+H"):
                error = precursor_error
                error_type = precursor_error_type
            charge = 1.0
            if ion_type.endswith("+2"):
                charge = 2.0
            if ion_type.endswith("+3"):
                charge = 3.0
            mz_ = mz_calc + (13.00335483778 - 12) / charge
            count = 0
            for ion in ions:
                if used_ions[count] == 0:
                    if error_type == "ppm":
                        if abs(ion[0] - mz_) < error * mz_ / 1e6:
                            mgf_pepmatch_ions_isotope1.append(
                                ["-" + ion_type, mz_, ion[0], ion[1]]
                            )
                            if ion_type.lstrip("-").startswith("b"):
                                sequence_evidence_b[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("b")
                                    )
                                    - 1
                                ] += ion[1]
                            if ion_type.lstrip("-").startswith("y"):
                                sequence_evidence_y[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("y")
                                    )
                                    - 1
                                ] += ion[1]
                            if used_ions[count] == 0 and frag_mz_min < ion[0]:
                                matched_intensity += ion[1]
                                if ion_type.lstrip("-").startswith("["):
                                    matched_intensity_parent += ion[1]
                                else:
                                    matched_intensity_not_parent += ion[1]
                            used_ions[count] = 1
                        elif mz_ < ion[0]:
                            break
                    else:
                        if abs(ion[0] - mz_) < error:
                            mgf_pepmatch_ions_isotope1.append(
                                ["-" + ion_type, mz_, ion[0], ion[1]]
                            )
                            if ion_type.lstrip("-").startswith("b"):
                                sequence_evidence_b[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("b")
                                    )
                                    - 1
                                ] += ion[1]
                            if ion_type.lstrip("-").startswith("y"):
                                sequence_evidence_y[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("y")
                                    )
                                    - 1
                                ] += ion[1]
                            if used_ions[count] == 0 and frag_mz_min < ion[0]:
                                matched_intensity += ion[1]
                                if ion_type.lstrip("-").startswith("["):
                                    matched_intensity_parent += ion[1]
                                else:
                                    matched_intensity_not_parent += ion[1]
                            used_ions[count] = 1
                        elif mz_ < ion[0]:
                            break
                count += 1
        mgf_pepmatch_ions_isotope2 = []
        for ion_type, mz_calc, mz_exp, intensity in mgf_pepmatch_ions_isotope1:
            error = fragment_error
            error_type = fragment_error_type
            if ion_type.lstrip("-").startswith("[M+H"):
                error = precursor_error
                error_type = precursor_error_type
            charge = 1.0
            if ion_type.endswith("+2"):
                charge = 2.0
            if ion_type.endswith("+3"):
                charge = 3.0
            mz_ = mz_calc + (13.00335483778 - 12) / charge
            count = 0
            for ion in ions:
                if used_ions[count] == 0:
                    if error_type == "ppm":
                        if abs(ion[0] - mz_) < error * mz_ / 1e6:
                            mgf_pepmatch_ions_isotope2.append(
                                ["-" + ion_type, mz_, ion[0], ion[1]]
                            )
                            if ion_type.lstrip("-").startswith("b"):
                                sequence_evidence_b[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("b")
                                    )
                                    - 1
                                ] += ion[1]
                            if ion_type.lstrip("-").startswith("y"):
                                sequence_evidence_y[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("y")
                                    )
                                    - 1
                                ] += ion[1]
                            if used_ions[count] == 0 and frag_mz_min < ion[0]:
                                matched_intensity += ion[1]
                                if ion_type.lstrip("-").startswith("["):
                                    matched_intensity_parent += ion[1]
                                else:
                                    matched_intensity_not_parent += ion[1]
                            used_ions[count] = 1
                        elif mz_ < ion[0]:
                            break
                    else:
                        if abs(ion[0] - mz_) < error:
                            mgf_pepmatch_ions_isotope2.append(
                                ["-" + ion_type, mz_, ion[0], ion[1]]
                            )
                            if ion_type.lstrip("-").startswith("b"):
                                sequence_evidence_b[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("b")
                                    )
                                    - 1
                                ] += ion[1]
                            if ion_type.lstrip("-").startswith("y"):
                                sequence_evidence_y[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("y")
                                    )
                                    - 1
                                ] += ion[1]
                            if used_ions[count] == 0 and frag_mz_min < ion[0]:
                                matched_intensity += ion[1]
                                if ion_type.lstrip("-").startswith("["):
                                    matched_intensity_parent += ion[1]
                                else:
                                    matched_intensity_not_parent += ion[1]
                            used_ions[count] = 1
                        elif mz_ < ion[0]:
                            break
                count += 1
        mgf_pepmatch_ions_isotope3 = []
        for ion_type, mz_calc, mz_exp, intensity in mgf_pepmatch_ions_isotope2:
            error = fragment_error
            error_type = fragment_error_type
            if ion_type.lstrip("-").startswith("[M+H"):
                error = precursor_error
                error_type = precursor_error_type
            charge = 1.0
            if ion_type.endswith("+2"):
                charge = 2.0
            if ion_type.endswith("+3"):
                charge = 3.0
            mz_ = mz_calc + (13.00335483778 - 12) / charge
            count = 0
            for ion in ions:
                if used_ions[count] == 0:
                    if error_type == "ppm":
                        if abs(ion[0] - mz_) < error * mz_ / 1e6:
                            mgf_pepmatch_ions_isotope3.append(
                                ["-" + ion_type, mz_, ion[0], ion[1]]
                            )
                            if ion_type.lstrip("-").startswith("b"):
                                sequence_evidence_b[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("b")
                                    )
                                    - 1
                                ] += ion[1]
                            if ion_type.lstrip("-").startswith("y"):
                                sequence_evidence_y[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("y")
                                    )
                                    - 1
                                ] += ion[1]
                            if used_ions[count] == 0 and frag_mz_min < ion[0]:
                                matched_intensity += ion[1]
                                if ion_type.lstrip("-").startswith("["):
                                    matched_intensity_parent += ion[1]
                                else:
                                    matched_intensity_not_parent += ion[1]
                            used_ions[count] = 1
                        elif mz_ < ion[0]:
                            break
                    else:
                        if abs(ion[0] - mz_) < error:
                            mgf_pepmatch_ions_isotope3.append(
                                ["-" + ion_type, mz_, ion[0], ion[1]]
                            )
                            if ion_type.lstrip("-").startswith("b"):
                                sequence_evidence_b[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("b")
                                    )
                                    - 1
                                ] += ion[1]
                            if ion_type.lstrip("-").startswith("y"):
                                sequence_evidence_y[
                                    int(
                                        ion_type.split("+")[0]
                                        .split("o")[0]
                                        .split("i")[0]
                                        .split("*")[0]
                                        .lstrip("-")
                                        .lstrip("y")
                                    )
                                    - 1
                                ] += ion[1]
                            if used_ions[count] == 0 and frag_mz_min < ion[0]:
                                matched_intensity += ion[1]
                                if ion_type.lstrip("-").startswith("["):
                                    matched_intensity_parent += ion[1]
                                else:
                                    matched_intensity_not_parent += ion[1]
                            used_ions[count] = 1
                        elif mz_ < ion[0]:
                            break
                count += 1
        for value in mgf_pepmatch_ions_isotope1:
            mgf_pepmatch_ions.append(value)
        for value in mgf_pepmatch_ions_isotope2:
            mgf_pepmatch_ions.append(value)
        for value in mgf_pepmatch_ions_isotope3:
            mgf_pepmatch_ions.append(value)
    return matched_intensity, matched_intensity_parent, matched_intensity_not_parent


def theoretical_matched_intensity(
    pep,
    scan_str,
    modifications,
    diff_dalton,
    charge,
    fragment,
    fragment_type,
    frag_mz_min,
    max_peaks_per_scan,
    mgf_table,
    current_target,
    masses,
):
    matched_intensity_max = 0
    mod = ""
    mod_with_maxmatchedintensity = str(modifications) + str(diff_dalton) + "@?"
    mod_with_maxmatchedintensity_u = ""

    mod_pos_mass = {}
    mod_pos = []
    if not (modifications == ""):
        for j in range(len(modifications[:-1].split(","))):
            mod_pos_mass[
                modifications.split(",")[j].split("@")[1][1:]
            ] = modifications.split(",")[j].split("@")[0]
            mod_pos.append(modifications.split(",")[j].split("@")[1][1:])
    sequence_evidence_b_with_maxmatchedintensity = np.zeros(len(pep))
    sequence_evidence_y_with_maxmatchedintensity = np.zeros(len(pep))

    for i in range(len(pep)):
        # i=10
        # make sure the dic coyied is not mutatable: https://www.peterbe.com/plog/be-careful-with-using-dict-to-create-a-copy
        mod_pos_mass_new = dict(mod_pos_mass)
        mod_pos_mass_unique = {}

        mgf_pepmatch_ions = []
        sequence_evidence_b = np.zeros(len(pep))
        sequence_evidence_y = np.zeros(len(pep))
        if str(i + 1) not in mod_pos:
            # print ('not'+str(i+1))
            mod_pos_mass_new[str(i + 1)] = str(diff_dalton)
            mod_pos_mass_unique[str(i + 1)] = str(diff_dalton)

        if str(i + 1) in mod_pos:
            # print ('yes'+str(i+1))
            mod_pos_mass_new[str(i + 1)] = str(
                float(mod_pos_mass[str(i + 1)]) + diff_dalton
            )
            mod_pos_mass_unique[str(i + 1)] = str(
                float(mod_pos_mass[str(i + 1)]) + diff_dalton
            )

        mod = ""
        mod_new_only = ""
        for j in mod_pos_mass_new:
            mod = mod + str(mod_pos_mass_new[j]) + "@" + str(pep[int(j) - 1]) + j + ","

        for u in mod_pos_mass_unique:
            mod_new_only = (
                mod_new_only
                + str(mod_pos_mass_unique[u])
                + "@"
                + str(pep[int(u) - 1])
                + u
                + ","
            )
            # print (i+1,mod)

        fragments = calc_peptide_fragments(fragments, pep, mod, charge, "by", masses)
        for c in range(charge, 0, -1):
            mass_ion_t = (
                (calc_peptide_mass(pep, mod) + c * masses["Proton"]) / (1.0 * c),
                "[M+H]+" + str(c),
            )
            fragments.append(mass_ion_t)
            # if '144.102' in mod:
            # 	mass_ion_t=(mass_new.calc_peptide_mass(pep, mod, masses)-144.10207+c*masses['Proton'])/(1.0*c),'[M+H-iTRAQ]+'+str(c)
            # 	fragments.append(mass_ion_t)
            for l in range(1, 3, 1):  # deamination and dehydration (1-2)
                l_text = ""
                if l > 1:
                    l_text = str(l) + "*"
                mass_ion_t = (
                    (
                        calc_peptide_mass(pep, mod)
                        - l * masses["H2O"]
                        + c * masses["Proton"]
                    )
                    / (1.0 * c),
                    "[M+H-" + l_text + "H2O]+" + str(c),
                )
                fragments.append(mass_ion_t)
                mass_ion_t = (
                    (
                        calc_peptide_mass(pep, mod)
                        - l * masses["NH3"]
                        + c * masses["Proton"]
                    )
                    / (1.0 * c),
                    "[M+H-" + l_text + "NH3]+" + str(c),
                )
                fragments.append(mass_ion_t)
        mz_error_fragment = fragment
        mz_error_fragment_type = fragment_type
        mz_error_precursor = 300.0
        mz_error_precursor_type = "ppm"

        if mz_error_fragment_type == "ppm":
            if mz_error_fragment == 60:
                mz_error_fragment = 100.0
                mz_error_precursor = 600.0
            else:
                mz_error_fragment = 20.0
        else:
            mz_error_fragment = 0.8
            mz_error_precursor = mz_error_fragment
            mz_error_precursor_type = mz_error_fragment_type

        (
            matched_intensity,
            matched_intensity_parent,
            matched_intensity_not_parent,
        ) = select_fragments_scan(
            scan_str,
            frag_mz_min,
            mgf_pepmatch_ions,
            sequence_evidence_b,
            sequence_evidence_y,
            fragments,
            mz_error_precursor,
            mz_error_precursor_type,
            mz_error_fragment,
            mz_error_fragment_type,
            max_peaks_per_scan,
            mgf_table,
            current_target,
        )
        # print i+1,matched_intensity
        if matched_intensity > matched_intensity_max:
            matched_intensity_max = matched_intensity
            mod_with_maxmatchedintensity = mod
            mod_with_maxmatchedintensity_u = mod_new_only
            # mgf_pepmatch_ions_with_maxmatchedintensity=mgf_pepmatch_ions
            sequence_evidence_b_with_maxmatchedintensity = sequence_evidence_b
            sequence_evidence_y_with_maxmatchedintensity = sequence_evidence_y
    # print  matched_intensity_max
    return (
        matched_intensity_max,
        mod_with_maxmatchedintensity,
        mod_with_maxmatchedintensity_u,
        sequence_evidence_b_with_maxmatchedintensity,
        sequence_evidence_y_with_maxmatchedintensity,
    )


@timeit
def main(args):
    """The main method for the agnostic search program

    Args:
        param args (param NamedTuple): A list of command line arguments supplied to the program

    Returns: 
        int: exit code of the software
    """

    # Steps:
    # 1. Parse the id file
    # 2. Load the MGF files
    # 3. Compare scans
    #   3.1 Calculate target scan information
    #   3.2 Calculate search window
    #   3.3 Compare within window
    # 4. Filter
    print(os.getcwd())
    print(args)

    print("LOADING QUERY MGF")
    mgf_queries = mgf_processor(os.path.abspath(os.path.expanduser(args.mgfQueryFile)))
    print("LOADING TARGET MGF")
    mgf_targets = mgf_processor(os.path.abspath(os.path.expanduser(args.mgfTargetFile)))
    print("LOADING ID MGF")
    mgf_filename = os.path.basename(args.mgfQueryFile).split(".")[0]
    df = pd.read_table(args.IdFile, index_col="scan")
    df = df.loc[df["Filename"] == mgf_filename]
    # FIXME: This replaces NA values with empty string (either make this clearer or remove if not necessary)
    # df = df.where((pd.notnull(df)), "")
    df = df.fillna("")
    df.index = df.index.astype(int)

    print(df, df.shape, df.index.dtype)
    q = list(filter(lambda query: query.scan in df.index, mgf_queries))

    data = []

    if not df.shape[0]:
        return -1
    # TODO: Flip target and query as there should be fewer queries than targets

    for target in mgf_targets:

        window = target.get_window()
        index_range = find_query_range(q, *window)

        best_match = 0.0
        bm = None
        for query in q[slice(*index_range)]:
            result = target.count_matches(query, args.mz_error, args.mz_error_type)
            # TODO: query matched threshold should be a cli argument
            if result.query_matched >= 0.3 and result.query_matched > best_match:
                best_match = result.query_matched
                bm = (query, result)

        if bm and bm[1].query_matched >= 0.3:
            query_scan_max = bm[0].scan
            peptide = df["peptide"].loc[query_scan_max]
            modifications = df["modifications"].loc[query_scan_max]
            expect = df["expect"].loc[query_scan_max]
            total_MS2_intensity = df["total_MS2_intensity"].loc[query_scan_max]
            proteins = df["proteins"].loc[query_scan_max]
            query_charge_max = int(df["charge"].loc[query_scan_max])
            query_scan_mz_max = bm[0].pepmass  # .m_zs[bm[0].top_idx]
            diff_dalton_max = bm[1].difference

            ####### add matched_theoretical_intensity ######
            frag_mz_min = 200
            # frag_mz_max = 0
            # frag_int_max = 0
            frag_int_sum = 1
            for idx in target.top_idx:
                if target.m_zs[idx] > frag_mz_min:
                    frag_int_sum += target.intensities[idx]
                    # scan_mgf_table = mgf_table_target[target_scan_string]
            (
                matched_intensity_max,
                mod_with_maxmatchedintensity,
                mod_with_maxmatchedintensity_u,
                sequence_evidence_b_with_maxmatchedintensity,
                sequence_evidence_y_with_maxmatchedintensity,
            ) = theoretical_matched_intensity(
                peptide,
                target.scan,
                modifications,
                diff_dalton_max,
                charge=target.charge,
                fragment=args.mz_error,
                fragment_type=args.mz_error_type,
                frag_mz_min=frag_mz_min,
                max_peaks_per_scan=args.max_peaks_per_target,
                mgf_table=mgf_targets,
                current_target=target,
                masses=MASSES,
            )
            largest_gap = calculate_largest_gap(
                peptide,
                sequence_evidence_b_with_maxmatchedintensity,
                sequence_evidence_y_with_maxmatchedintensity,
            )
            largest_gap_per = largest_gap * 1.0 / len(peptide)
            ####### add matched_theoretical_intensity ######
            data.append(
                [
                    str(os.path.basename(args.mgfQueryFile)),
                    # str(mgf_query_file_only),
                    # str(mgf_target_file_only),
                    str(os.path.basename(args.mgfTargetFile)),
                    str("{:.4f}".format(diff_dalton_max)),
                    # str("{:.1f}".format(float(diff_int_max))),
                    str(
                        "{:.1f}".format(
                            float(
                                math.copysign(
                                    math.ceil(abs(diff_dalton_max)), diff_dalton_max
                                )
                            )
                        )
                    ),
                    bm[0].scan,
                    str(query_scan_mz_max),
                    str(query_charge_max),
                    target.scan,
                    # str(target_scan_mz),
                    str(target.pepmass),
                    # str(mgf_charge_target[target_scan_string]),
                    str(target.charge),
                    # str(matches_max),
                    str(bm[1].count),
                    # str(matched_query_max),
                    str(bm[1].query_matched),
                    # str(matched_query_max * matched_target_max),
                    str(bm[1].query_matched * bm[1].target_matched),
                    # str(sum_query_max + sum_target_max),
                    str(bm[1].query_log_sum + bm[1].target_log_sum),
                    str(peptide),
                    str(modifications),
                    str(expect),
                    str(total_MS2_intensity),
                    str(proteins),
                    "{:.4f}".format(matched_intensity_max / (1.0 * frag_int_sum)),
                    str(largest_gap),
                    str(largest_gap_per),
                    str(mod_with_maxmatchedintensity),
                    str(mod_with_maxmatchedintensity_u),
                ]
            )

    out_df = pd.DataFrame(
        data,
        columns=[
            "MGF_query_file",
            "MGF_target_file",
            "Diff_dalton",
            "Diff_dalton_bin",
            "Query_scan",
            "Query_scan_mz",
            "Query_scan_charge",
            "Target_scan",
            "Target_scan_mz",
            "Target_scan_charge",
            "Matches",
            "Matched_query",
            "Matched_intensity_product",
            "Sum_log_intensity",
            "Peptide",
            "Modifications",
            "Expect",
            "Total_Query_MS2_intensity",
            "Proteins",
            "Matched_peptide_intensity_max",
            "Largest_gap",
            "Largest_gap_percent",
            "Full_mod",
            "Unique_mod",
        ],
    )
    print(out_df)
    out_df.to_csv("test.csv", index=False)
    # results = pd.DataFrame(
    #     r,
    #     columns=[
    #         "count",
    #         "target_matched",
    #         "query_matched",
    #         "target_log",
    #         "query_log",
    #     ],
    # )

    # FIXME: abs difference >= 10 : This should be a cli option

    return 0
