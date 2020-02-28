import os
import pandas as pd

from src.agnostic_search.mgf.processer import process as mgf_processor


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

    mgf_queries = mgf_processor(os.path.abspath(args.mgfQueryFile))
    mgf_targets = mgf_processor(os.path.abspath(args.mgfTargetFile))
    mgf_filename = os.path.basename(args.mgfQueryFile).split(".")[0]
    df = pd.read_table(args.IdFile, index_col="scan")
    df = df.loc[df["Filename"] == mgf_filename]
    # FIXME: This replaces NA values with empty string (either make this clearer or remove if not necessary)
    df = df.where((pd.notnull(df)), "")
    df = df.fillna("")

    print(df, df.shape)
    count = 0
    q = list(filter(lambda x: x.scan in df.index, mgf_queries))
    if df.shape[0]:
        for target in mgf_targets:

            window = target.get_window()
            index_range = find_query_range(q, *window)
            for query in q[slice(*index_range)]:
                if query.scan in df.index:
                    count += 1
                    # print("-----------")
                    # print("target")
                    # print(target)
                    # print("query")
                    # print(query)
    print(count)
    return 0
