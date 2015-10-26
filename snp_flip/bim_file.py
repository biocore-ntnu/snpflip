from itertools import groupby

import pandas as pd


def read_bim_file(bim_file):
    bim_table = pd.read_table(bim_file, sep=r"\s+", header=None,
                              names=["chromosome", "snp_name",
                                     "genetic_distance", "position",
                                     "allele_1", "allele_2"],
                              dtype={"chromosome": str})

    bim_table["position"] = bim_table.position - 1

    return bim_table.set_index(["chromosome", "position"])
