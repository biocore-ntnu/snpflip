import pandas as pd


def read_bim_file(bim_file):
    return pd.read_table(bim_file, sep=r"\s+", header=None,
                         names=["chromosome", "snp_name", "genetic_distance",
                                "position", "allele_1", "allele_2"])
