from __future__ import print_function

import re
from sys import stderr

from operator import itemgetter

import pandas as pd

from pyfaidx import Fasta



def get_reference_genome_data(bim_table, genome_fasta):

    pyfaidx_genome = Fasta(genome_fasta, as_raw=True, read_ahead=100000)
    genome_dict = convert_fa_chromosome_names(pyfaidx_genome)

    chromosome_data = []
    for chromosome, nucleotides in genome_dict.items():

        positions = bim_table[bim_table["chromosome"] == chromosome]["position"]
        _get_snps = itemgetter(*list(positions))

        snp_nucleotides = [snp.upper() for snp in _get_snps(str(nucleotides))]
        chromosome_data.append(pd.Series(snp_nucleotides))

        _check_for_N(snp_nucleotides, chromosome)

    return pd.concat(chromosome_data)

def _check_for_N(snps, chromosome):

    nb_ns = sum([1 for snp in snps if snp == "N"])
    if nb_ns != 0:
        print("There were {nb_ns} 'N' nucleotides in chromosome"
              "{chromosome}.".format(**vars()), file=stderr)


def convert_fa_chromosome_names(fasta_dict):

    return {_convert_fa_name(k): fasta_dict[k] for k in fasta_dict.keys()
            if _is_valid_chromosome(k)}


def _is_valid_chromosome(name):

    for c in "_.":
        if c in name.split(None, 1)[0]:
            return False

    return True


def _convert_fa_name(name):

    bim_name = re.compile(r"(chr)?(?P<chr_symbol>\d+|X|Y|M)(?:T)?")
    result = bim_name.match(name).group("chr_symbol")

    return result
