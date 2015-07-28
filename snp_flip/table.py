import pandas as pd

from snp_flip.reference_genome import get_reference_genome_data

from snp_flip.bim_file import read_bim_file

def create_snp_table(bim_file, fa_file):

    bim_table = read_bim_file(bim_file)

    reference_genome_data = get_reference_genome_data(bim_table, fa_file)

    snp_table = pd.concat([bim_table, reference_genome_data], axis=1)

    return snp_table
