import math
import numpy as np
import os
import pandas as pd
import pprint
import time

from src.sampei.mgf.processer import process as mgf_processor
from src.sampei.masses.masses import (
    MASSES,
    calc_peptide_mass,
    calc_peptide_fragments,
)
from src.sampei.masses.calculations import select_fragments_scan


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print("%r  %2.2f s" % (method.__name__, (te - ts)))
        return result

    return timed


def find_range(mgf_list, low, high):
    lower_bound = 0
    upper_bound = len(mgf_list)
    while lower_bound < upper_bound:
        mid = lower_bound + (upper_bound - lower_bound) // 2
        if low <= mgf_list[mid].pepmass * mgf_list[mid].charge:
            upper_bound = mid
        else:
            lower_bound = mid + 1
    low_idx = lower_bound
    lower_bound = 0
    upper_bound = len(mgf_list)

    while lower_bound < upper_bound:
        mid = lower_bound + (upper_bound - lower_bound) // 2
        if high < mgf_list[mid].pepmass * mgf_list[mid].charge:
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


def theoretical_matched_intensity(
    pep,
    modifications,
    diff_dalton,
    charge,
    fragment,
    fragment_type,
    frag_mz_min,
    max_peaks_per_scan,
    current_target,
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

        fragments = calc_peptide_fragments(pep, mod, charge, "by")
        for c in range(charge, 0, -1):
            mass_ion_t = (
                (calc_peptide_mass(pep, mod) + c * MASSES["Proton"]) / (1.0 * c),
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
                        - l * MASSES["H2O"]
                        + c * MASSES["Proton"]
                    )
                    / (1.0 * c),
                    "[M+H-" + l_text + "H2O]+" + str(c),
                )
                fragments.append(mass_ion_t)
                mass_ion_t = (
                    (
                        calc_peptide_mass(pep, mod)
                        - l * MASSES["NH3"]
                        + c * MASSES["Proton"]
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
        ions = list(
            zip(
                current_target.m_zs[current_target.top_idx],
                current_target.intensities[current_target.top_idx],
            )
        )
        (
            matched_intensity,
            sequence_evidence_b,
            sequence_evidence_y,
            mgf_pepmatch_ions,
        ) = select_fragments_scan(
            frag_mz_min,
            fragments,
            mz_error_precursor,
            mz_error_precursor_type,
            mz_error_fragment,
            mz_error_fragment_type,
            max_peaks_per_scan,
            ions,
            len(pep),
        )
        # print i+1,matched_intensity
        if matched_intensity > matched_intensity_max:
            matched_intensity_max = matched_intensity
            mod_with_maxmatchedintensity = mod
            mod_with_maxmatchedintensity_u = mod_new_only
            mgf_pepmatch_ions_with_maxmatchedintensity = mgf_pepmatch_ions
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
    # TODO: Add a validation step
    print(os.getcwd())
    print(args)

    print("LOADING QUERY MGF")
    mgf_queries = mgf_processor(
        os.path.abspath(os.path.expanduser(args.mgf_query_file)),
        args.max_peaks_per_scan,
    )

    print("LOADING TARGET MGF")
    mgf_targets = mgf_processor(
        os.path.abspath(os.path.expanduser(args.mgf_target_file)),
        args.max_peaks_per_scan,
    )
    print("LOADING ID MGF")
    mgf_filename = os.path.basename(args.mgf_query_file).split(".")[0]
    df = pd.read_table(args.id_file, index_col="scan")
    df = df.loc[df["Filename"] == mgf_filename]
    df = df.fillna("")
    df.index = df.index.astype(int)

    print(df, df.shape, df.index.dtype)
    filtered_queries = list(filter(lambda query: query.scan in df.index, mgf_queries))

    data = []

    if not df.shape[0]:
        return -1
    for t_idx, target in enumerate(mgf_targets):

        window = target.get_window()
        index_range = find_range(filtered_queries, *window)

        best_match = 0.0
        bm = None
        if t_idx % 1000 == 0:
            print(t_idx)
        if index_range[0] == index_range[1]:
            continue

        for query in filtered_queries[slice(*index_range)]:
            result = target.count_matches(
                query, args.fragment_mass_error, args.error_type
            )
            if (
                result.query_matched >= args.matched_query_intensity
                and result.query_matched > best_match
            ):
                best_match = result.query_matched
                bm = (query, result)

        if bm and bm[1].query_matched >= args.matched_query_intensity:
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
                modifications,
                diff_dalton_max,
                charge=target.charge,
                fragment=args.fragment_mass_error,
                fragment_type=args.error_type,
                frag_mz_min=frag_mz_min,
                max_peaks_per_scan=args.max_peaks_per_scan,
                current_target=target,
            )
            largest_gap = calculate_largest_gap(
                peptide,
                sequence_evidence_b_with_maxmatchedintensity,
                sequence_evidence_y_with_maxmatchedintensity,
            )
            largest_gap_per = largest_gap / len(peptide)
            ####### add matched_theoretical_intensity ######
            data.append(
                [
                    str(os.path.basename(args.mgf_query_file)),
                    # str(mgf_query_file_only),
                    # str(mgf_target_file_only),
                    str(os.path.basename(args.mgf_target_file)),
                    str("{:.4f}".format(diff_dalton_max)),
                    # str("{:.1f}".format(float(diff_int_max))),
                    float(
                        math.copysign(math.ceil(abs(diff_dalton_max)), diff_dalton_max)
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
                    matched_intensity_max / frag_int_sum,
                    largest_gap,
                    largest_gap_per,
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

    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)
    output_file = "{}/query-{}.target-{}.error_unit-{}.error-{}.peaks-{}".format(
        args.output_directory,
        os.path.basename(args.mgf_query_file[:-4]),
        os.path.basename(args.mgf_target_file[:-4]),
        args.error_type,
        args.fragment_mass_error,
        args.max_peaks_per_scan,
    )
    if not args.no_filter:
        output_file += ".lgp-{}.mpi-{}".format(
            args.largest_gap_percent, args.min_diff_dalton_bin
        )

        print("Filtering unwanted mass shift/modifications...")
        #     ###remove unexpected modifications produced by x!tandem####
        unexpected_modifications = (
            out_df["Modifications"].str.contains("-17.02655")
            | out_df["Modifications"].str.contains("42.01057")
            | out_df["Modifications"].str.contains("-18.01056")
        ).values.sum()

        out_df = out_df[
            ~(
                out_df["Modifications"].str.contains("-17.02655")
                | out_df["Modifications"].str.contains("42.01057")
                | out_df["Modifications"].str.contains("-18.01056")
            )
        ]
        num_less_gap_percent = (
            (out_df["Largest_gap_percent"] <= args.largest_gap_percent).values.sum()
            & (
                out_df["Matched_peptide_intensity_max"]
                >= args.matched_peptide_intensity
            )
        ).values.sum()
        out_df_filtered = out_df[
            (out_df["Largest_gap_percent"] <= args.largest_gap_percent)
            & (
                out_df["Matched_peptide_intensity_max"]
                >= args.matched_peptide_intensity
            )
        ]
        out_df = out_df_filtered[
            abs(out_df_filtered["Diff_dalton_bin"]) > args.min_diff_dalton_bin
        ]
        if args.write_intermediate:
            out_df_filtered.to_csv(
                output_file + ".tab", sep="\t", index=False,
            )
            out_df.to_csv(
                output_file
                + ".dalton_bin-{}".format(args.min_diff_dalton_bin)
                + ".tab",
                sep="\t",
                index=False,
            )
        print("Done filtering unwanted mass shift/modifications...")

        if args.xtandem_xml:
            output_file += ".noxtandem"
            print("Reading xtandem ouput..")
            xtandem_all = pd.read_table(args.xtandem_xml, sep="\t")
            out_df = out_df[
                ~out_df.set_index(["Target_scan"]).index.isin(
                    xtandem_all.set_index(["scan"]).index
                )
            ]
    out_df.to_csv(output_file + ".tab", sep="\t", index=False)
    print(out_df)
    return 0
