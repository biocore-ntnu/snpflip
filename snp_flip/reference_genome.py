from __future__ import print_function

import re
from sys import stderr

from operator import itemgetter

from pyfaidx import Fasta

def get_reference_genome_data(snp_positions, chromosome, genome_fasta):

    _get_snps = itemgetter(*snp_positions)
    chrom_data = Fasta(genome_fasta, as_raw=True, read_ahead=100000)[chromosome]
    snp_nucleotides = [snp.upper() for snp in _get_snps(chrom_data)]

    _check_for_N(snp_nucleotides, chromosome)

    return snp_nucleotides

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
