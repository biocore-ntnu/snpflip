from itertools import groupby

import pandas as pd


def read_bim_file(bim_file):
    bim_table = pd.read_table(bim_file, sep=r"\s+", header=None,
                              names=["chromosome", "snp_name",
                                     "genetic_distance", "position",
                                     "allele_1", "allele_2"],
                              dtype={"chromosome": str})

    _check_if_rows_grouped_by_chromosomes(bim_table["chromosome"])

    return bim_table


def _check_if_rows_grouped_by_chromosomes(chromosomes):


    nb_chromosomes = len(chromosomes.drop_duplicates())
    uniq_chromosomes = len(list(groupby(chromosomes, lambda x: x)))

    if nb_chromosomes != uniq_chromosomes:
        raise ValueError("The rows in the bim file are not grouped by " \
                         "chromosomes. Exiting.")
