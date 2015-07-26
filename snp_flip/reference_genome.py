from __future__ import print_function

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
