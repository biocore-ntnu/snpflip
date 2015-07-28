import pandas as pd

from snp_flip.reference_genome import get_reference_genome_data

from snp_flip.bim_file import read_bim_file

def create_snp_table(bim_file, fa_file):

    bim_table = read_bim_file(bim_file)

    reference_genome_data = get_reference_genome_data(bim_table, fa_file)

    snp_table = pd.concat([bim_table, reference_genome_data], axis=1)

    return snp_table


def find_strand(bim_table):

    forward_mask = (bim_table["reference"] == bim_table["allele_1"]) | \
                   (bim_table["reference"] == bim_table["allele_2"])

    reverse_mask = (bim_table["reference_rev"] == bim_table["allele_1"]) | \
                   (bim_table["reference_rev"] == bim_table["allele_2"])

    ambiguous_mask = (forward_mask & reverse_mask) | \
                     (~(forward_mask | reverse_mask))

    new_series = bim_table["reference"].copy()
    new_series[forward_mask] = "forward"
    new_series[reverse_mask] = "reverse"
    new_series[ambiguous_mask] = "ambiguous"
    new_series.name = "strand"

    return pd.concat([bim_table, new_series], axis=1)
