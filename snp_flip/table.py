import pandas as pd

from snp_flip.reference_genome import get_reference_genome_data

from snp_flip.bim_file import read_bim_file

def create_snp_table(bim_file, fa_file):

    bim_table = read_bim_file(bim_file)

    reference_genome_data = get_reference_genome_data(bim_table, fa_file)

    snp_table = bim_table.join(reference_genome_data)

    return snp_table


def find_strand(bim_table):

    forward_mask = (bim_table["reference"] == bim_table["allele_1"]) | \
                   (bim_table["reference"] == bim_table["allele_2"])

    reverse_mask = (bim_table["reference_rev"] == bim_table["allele_1"]) | \
                   (bim_table["reference_rev"] == bim_table["allele_2"])

    N_mask = bim_table["reference"] == "N"

    ambiguous_mask = (forward_mask & reverse_mask) | \
                     (~(forward_mask | reverse_mask)) | \
                     N_mask


    strand_series = bim_table["reference"].copy()
    strand_series[forward_mask] = "forward"
    strand_series[reverse_mask] = "reverse"
    strand_series[ambiguous_mask] = "ambiguous"
    strand_series.name = "strand"

    return pd.concat([bim_table, strand_series], axis=1)
