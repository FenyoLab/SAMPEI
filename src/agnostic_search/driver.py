import os

from src.agnostic_search.mgf_processer import main as mgf


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
    mgf(os.path.abspath(args.mgfQueryFile))

    print(args)
    return 0
