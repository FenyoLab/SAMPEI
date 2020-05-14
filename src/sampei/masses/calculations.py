import numpy as np

ION_MATCHING = {
    "ppm": lambda ion, frag, error: abs(ion - frag) < error * frag / 1e6,
    "dalton": lambda ion, frag, error: abs(ion - frag) < error,
}


def select_fragments_scan(
    frag_mz_min,
    fragments,
    precursor_error,
    precursor_error_type,
    fragment_error,
    fragment_error_type,
    num_intense,
    ions,
    pep_len,
):
    sequence_evidence = {"b": np.zeros(pep_len), "y": np.zeros(pep_len)}
    mgf_pepmatch_ions = []
    matched_intensity = 0.0
    matched_intensity_not_parent = 0.0
    matched_intensity_parent = 0.0
    # ions = mgf_table[scan]

    # sort ions for getting most intense
    if num_intense > 0:
        ions.sort(
            key=lambda tup: tup[1], reverse=True
        )  # sorts descending, by intensity
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
            if ION_MATCHING[error_type](ion[0], frag_mz[0], error):
                mgf_pepmatch_ions.append([frag_mz[1], frag_mz[0], ion[0], ion[1]])
                b_y_ion = frag_mz[1].lstrip("-")[0]
                if b_y_ion in "by":
                    sequence_evidence[b_y_ion][
                        int(
                            frag_mz[1]
                            .split("+")[0]
                            .split("o")[0]
                            .split("i")[0]
                            .split("*")[0]
                            .lstrip("-")
                            .lstrip(b_y_ion)
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
    if fragment_error_type == "ppm":
        mgf_pepmatch_ions_isotope1 = get_pepmatch_isotope(
            mgf_pepmatch_ions,
            fragment_error,
            fragment_error_type,
            precursor_error,
            precursor_error_type,
            ions,
            used_ions,
            sequence_evidence,
            frag_mz_min,
            matched_intensity,
            matched_intensity_parent,
            matched_intensity_not_parent,
        )
        mgf_pepmatch_ions_isotope2 = get_pepmatch_isotope(
            mgf_pepmatch_ions_isotope1,
            fragment_error,
            fragment_error_type,
            precursor_error,
            precursor_error_type,
            ions,
            used_ions,
            sequence_evidence,
            frag_mz_min,
            matched_intensity,
            matched_intensity_parent,
            matched_intensity_not_parent,
        )
        mgf_pepmatch_ions_isotope3 = get_pepmatch_isotope(
            mgf_pepmatch_ions_isotope2,
            fragment_error,
            fragment_error_type,
            precursor_error,
            precursor_error_type,
            ions,
            used_ions,
            sequence_evidence,
            frag_mz_min,
            matched_intensity,
            matched_intensity_parent,
            matched_intensity_not_parent,
        )
        for value in mgf_pepmatch_ions_isotope1:
            mgf_pepmatch_ions.append(value)
        for value in mgf_pepmatch_ions_isotope2:
            mgf_pepmatch_ions.append(value)
        for value in mgf_pepmatch_ions_isotope3:
            mgf_pepmatch_ions.append(value)
    return (
        matched_intensity,
        # matched_intensity_parent,
        # matched_intensity_not_parent,
        sequence_evidence["b"],
        sequence_evidence["y"],
        mgf_pepmatch_ions,
    )


def get_pepmatch_isotope(
    mgf_pepmatch_ions,
    fragment_error,
    fragment_error_type,
    precursor_error,
    precursor_error_type,
    ions,
    used_ions,
    sequence_evidence,
    frag_mz_min,
    matched_intensity,
    matched_intensity_parent,
    matched_intensity_not_parent,
):
    mgf_pepmatch_ions_isotope = []
    for ion_type, mz_calc, _, _ in mgf_pepmatch_ions:
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
                if ION_MATCHING[error_type](ion[0], mz_, error):
                    mgf_pepmatch_ions_isotope.append(
                        ["-" + ion_type, mz_, ion[0], ion[1]]
                    )
                    b_y_ion = ion_type.lstrip("-")[0]
                    if b_y_ion in "by":
                        sequence_evidence[b_y_ion][
                            int(
                                ion_type.split("+")[0]
                                .split("o")[0]
                                .split("i")[0]
                                .split("*")[0]
                                .lstrip("-")
                                .lstrip(b_y_ion)
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
    return mgf_pepmatch_ions_isotope
