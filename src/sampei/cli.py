import argparse

from src.sampei.driver import main as driver


def main():
    parser = argparse.ArgumentParser(
        description="SAMPEI is a searching method leveraging high quality query spectra within the same or different dataset to assign target spectra with peptide sequence and undefined modification (mass shift).",
        prog='SAMPEI'
    
    )
    parser.add_argument(
        "mgf_query_file",
        type=str,
        help="Query mgf file with full path containing query scans have been identified by DB search",
    )
    parser.add_argument(
        "mgf_target_file",
        type=str,
        help="Target mgf file with full path containing target scans with undefined modifications",
    )
    parser.add_argument(
        "id_file",
        type=str,
        help="File in which query scans have been identified by DB search",
    )
    parser.add_argument(
        "--max-peaks-per-scan", type=int, default=20, help="(default = %(default)s)"
    )
    # parser.add_argument("--precursor-mass-error", type=int, default=20, help="")
    parser.add_argument(
        "--fragment-mass-error", type=int, default=20, help="(default = %(default)s)"
    )
    parser.add_argument(
        "--matched-query-intensity",
        type=float,
        default=0.3,
        help="The percentage of MS2 intensity of query scan matched to target scan over the summation of total MS2 intensity in the query scan. (default = %(default)s)",
    )
    parser.add_argument(
        "--error-type",
        type=str,
        choices=["ppm", "dalton"],
        default="ppm",
        help="(default = %(default)s)",
    )
    parser.add_argument(
        "-O",
        "--output-directory",
        type=str,
        default="output",
        help="Full path to the directory where output is stored. If this directory does not exist it will be created. (default = %(default)s)",
    )
    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Disable the filter and keep DB identified scans in the target mgf file",
    )

    parser.add_argument(
        "--largest-gap-percent",
        type=float,
        default=0.4,
        help="The percentage of the largest consecutive b/y ion missing over the length of the peptide sequence. (default = %(default)s)",
    )
    parser.add_argument(
        "--matched-peptide-intensity",
        type=float,
        default=0.5,
        help="The percentage of MS2 intensity of target scan matched to the theoretical fragments of peptide sequence over the summation of total MS2 intensity in the target scan. (default = %(default)s)",
    )
    parser.add_argument(
        "--min-diff-dalton-bin",
        type=int,
        default=10,
        help="The absolute minimum dalton difference between the query scan and the target scan. (default = %(default)s)",
    )
    parser.add_argument(
        "--xtandem-xml",
        type=str,
        help="The path to an X!tandem xml file which will be used to filter the results.",
    )
    parser.add_argument(
        "--write-intermediate",
        action="store_true",
        help="Write files for each step of filtering.",
    )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.9')

    args = parser.parse_args()
    driver(args)


if __name__ == "__main__":
    main()
