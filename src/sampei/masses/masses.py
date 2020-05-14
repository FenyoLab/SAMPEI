import re
from string import ascii_letters

MASSES = {
    "Proton": 1.007276,
    "H": 1.007825035,
    "O": 15.99491463,
    "N": 14.003074,
    "C": 12.0,
    "C13": 13.00335483778,
    "S": 31.9720707,
    "P": 30.973762,
    "aa:J": 0.0,
    "aa:O": 0.0,
    "aa:U": 0.0,
    "aa:X": 0.0,
}


def calc_mass(composition):
    mass = 0.0
    formula_p = re.compile("([A-Z][a-z]?)([0-9]*)")
    comp = formula_p.findall(composition)
    for match in comp:
        atom, number = match
        if len(number) == 0:
            number = 1
        mass += int(number) * MASSES[atom]
    return mass


MASSES["H2O"] = 2 * MASSES["H"] + MASSES["O"]
MASSES["NH3"] = MASSES["N"] + 3 * MASSES["H"]
MASSES["HPO3"] = MASSES["H"] + MASSES["P"] + 3 * MASSES["O"]
MASSES["H3PO4"] = 3 * MASSES["H"] + MASSES["P"] + 4 * MASSES["O"]
MASSES["b"] = MASSES["Proton"]
MASSES["c"] = MASSES["N"] + 3 * MASSES["H"] + MASSES["Proton"]
MASSES["y"] = MASSES["O"] + 2 * MASSES["H"] + MASSES["Proton"]
MASSES["z"] = (
    2 * MASSES["H"] + MASSES["O"] - 2 * MASSES["H"] - MASSES["N"] + MASSES["Proton"]
)
MASSES["aa:A"] = calc_mass("C3H5ON")
MASSES["aa:B"] = calc_mass("C4H6O2N2")  # samee as N
MASSES["aa:C"] = calc_mass("C3H5ONS")
MASSES["aa:D"] = calc_mass("C4H5O3N")
MASSES["aa:E"] = calc_mass("C5H7O3N")
MASSES["aa:F"] = calc_mass("C9H9ON")
MASSES["aa:G"] = calc_mass("C2H3ON")
MASSES["aa:H"] = calc_mass("C6H7ON3")
MASSES["aa:I"] = calc_mass("C6H11ON")
MASSES["aa:K"] = calc_mass("C6H12ON2")
MASSES["aa:L"] = calc_mass("C6H11ON")
MASSES["aa:M"] = calc_mass("C5H9ONS")
MASSES["aa:N"] = calc_mass("C4H6O2N2")
MASSES["aa:P"] = calc_mass("C5H7ON")
MASSES["aa:Q"] = calc_mass("C5H8O2N2")
MASSES["aa:R"] = calc_mass("C6H12ON4")
MASSES["aa:S"] = calc_mass("C3H5O2N")
MASSES["aa:T"] = calc_mass("C4H7O2N")
MASSES["aa:V"] = calc_mass("C5H9ON")
MASSES["aa:W"] = calc_mass("C11H10ON2")
MASSES["aa:Y"] = calc_mass("C9H9O2N")
MASSES["aa:Z"] = calc_mass("C5H8O2N2")  # Same as Q


def calc_peptide_mass(peptide, modifications):
    mass = 0.0
    mod_dict = {}
    if len(modifications) > 0:
        mods = modifications.split(",")
        for m in mods:
            if "@" in m:
                mod, aa = m.split("@")
                if (
                    mod[0] in ascii_letters
                ):  # fix this! - to account for possible +/- before the formula
                    mod = calc_mass(mod)
                if len(aa) > 1:
                    if aa[0] in ascii_letters and aa[1] in [
                        "0",
                        "1",
                        "2",
                        "3",
                        "4",
                        "5",
                        "6",
                        "7",
                        "8",
                        "9",
                    ]:
                        aa = aa[1:]
                mod_dict[aa] = float(mod)
    for i, p in enumerate(peptide):
        mass += MASSES["aa:" + p]
        if str(i + 1) in mod_dict:
            mass += mod_dict[str(i + 1)]
        if p in mod_dict:
            mass += mod_dict[p]
    mass += MASSES["H2O"]
    return mass


def calc_peptide_fragments(peptide, modifications, charge, ion_types):
    fragments = []
    mass = 0.0
    mod_dict = {}
    if len(modifications) > 0:
        mods = modifications.split(",")
        for m in mods:
            if "@" in m:
                mod, aa = m.split("@")
                if (
                    mod[0] in ascii_letters
                ):  # fix this! - to account for possible +/- before the formula
                    mod = calc_mass(mod)
                if len(aa) > 1:
                    if aa[0] in ascii_letters and aa[1] in [
                        "0",
                        "1",
                        "2",
                        "3",
                        "4",
                        "5",
                        "6",
                        "7",
                        "8",
                        "9",
                    ]:
                        aa = aa[1:]
                mod_dict[aa] = float(mod)
    for ion_type in ion_types:
        mass = MASSES[ion_type]
        if ion_type in "bc":
            for i, p in enumerate(peptide):
                mass += MASSES["aa:" + p]
                if str(i + 1) in mod_dict:
                    mass += mod_dict[str(i + 1)]
                if p in mod_dict:
                    mass += mod_dict[p]
                if 0 <= i < len(peptide) - 1:
                    for k in range(1):  # (3):
                        prefix = ""
                        if k > 0:
                            prefix = "-"
                        mass_ = mass + k * (MASSES["C13"] - MASSES["C"])
                        mass_ion_t = mass_, prefix + ion_type + str(i + 1)
                        fragments.append(mass_ion_t)
                        if "1" in mod_dict:
                            if abs(mod_dict["1"] - 144.10207) < 0.01:
                                mass_ion_t = (
                                    mass_ - 144.10207,
                                    prefix + ion_type + str(i + 1) + "i",
                                )
                                fragments.append(mass_ion_t)
                        for l in range(1, 2, 1):
                            mass_ion_t = (
                                mass_ - l * MASSES["H2O"],
                                prefix + ion_type + str(i + 1) + "o" * l,
                            )
                            fragments.append(mass_ion_t)
                            mass_ion_t = (
                                mass_ - l * MASSES["NH3"],
                                prefix + ion_type + str(i + 1) + "*" * l,
                            )
                            fragments.append(mass_ion_t)
                        if int(charge) > 2:
                            mass_ion_t = (
                                (mass_ + MASSES["Proton"]) / 2.0,
                                prefix + ion_type + str(i + 1) + "+2",
                            )
                            fragments.append(mass_ion_t)
                            for l in range(1, 2, 1):
                                mass_ion_t = (
                                    (mass_ - l * MASSES["H2O"] + MASSES["Proton"])
                                    / 2.0,
                                    prefix + ion_type + str(i + 1) + "o" * l + "+2",
                                )
                                fragments.append(mass_ion_t)
                                mass_ion_t = (
                                    (mass_ - l * MASSES["NH3"] + MASSES["Proton"])
                                    / 2.0,
                                    prefix + ion_type + str(i + 1) + "*" * l + "+2",
                                )
                                fragments.append(mass_ion_t)
        else:
            for i in range(len(peptide) - 1, -1, -1):
                mass += MASSES["aa:" + peptide[i]]
                if str(i + 1) in mod_dict:
                    mass += mod_dict[str(i + 1)]
                if peptide[i] in mod_dict:
                    mass += mod_dict[peptide[i]]
                if 0 < i < len(peptide):
                    for k in range(1):  # (3):
                        prefix = ""
                        if k > 0:
                            prefix = "-"
                        mass_ = mass + k * (MASSES["C13"] - MASSES["C"])
                        mass_ion_t = mass_, prefix + ion_type + str(len(peptide) - i)
                        fragments.append(mass_ion_t)
                        if str(len(peptide)) in mod_dict:
                            if abs(mod_dict[str(len(peptide))] - 144.10207) < 0.01:
                                mass_ion_t = (
                                    mass_ - 144.10207,
                                    prefix + ion_type + str(len(peptide) - i) + "i",
                                )
                                fragments.append(mass_ion_t)
                        for l in range(1, 2, 1):
                            mass_ion_t = (
                                mass_ - l * MASSES["H2O"],
                                prefix + ion_type + str(len(peptide) - i) + "o" * l,
                            )
                            fragments.append(mass_ion_t)
                            mass_ion_t = (
                                mass_ - l * MASSES["NH3"],
                                prefix + ion_type + str(len(peptide) - i) + "*" * l,
                            )
                            fragments.append(mass_ion_t)
                        if int(charge) > 2:
                            mass_ion_t = (
                                (mass_ + MASSES["Proton"]) / 2.0,
                                prefix + ion_type + str(len(peptide) - i) + "+2",
                            )
                            fragments.append(mass_ion_t)
                            for l in range(1, 2, 1):
                                mass_ion_t = (
                                    (mass_ - l * MASSES["H2O"] + MASSES["Proton"])
                                    / 2.0,
                                    prefix
                                    + ion_type
                                    + str(len(peptide) - i)
                                    + "o" * l
                                    + "+2",
                                )
                                fragments.append(mass_ion_t)
                                mass_ion_t = (
                                    (mass_ - l * MASSES["NH3"] + MASSES["Proton"])
                                    / 2.0,
                                    prefix
                                    + ion_type
                                    + str(len(peptide) - i)
                                    + "*" * l
                                    + "+2",
                                )
                                fragments.append(mass_ion_t)
    return fragments
