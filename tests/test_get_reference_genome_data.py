import pytest

from numpy import array_equal
import pandas as pd


from ebs.imports import StringIO


from snp_flip.reference_genome import get_reference_genome_data

def describe_get_reference_genome_data():

    def test_get_data(bim_df, genome_fasta, expected_result):
        snp_nucleotides = get_reference_genome_data(bim_df, str(genome_fasta))
        print(snp_nucleotides)
        print(type(snp_nucleotides))
        print(expected_result)
        print(type(expected_result))
        assert array_equal(expected_result, snp_nucleotides)


@pytest.fixture
def expected_result():
    return pd.read_table(StringIO("""A  T\nT   A\nC   G\nA   T\n"""),
                         header=None, sep=r"\s+")


@pytest.fixture
def genome_fasta(tmpdir):

    genome_fa = tmpdir.join("genome.fa")
    genome_fa.write(u""">chr1\nACGTTCCG\n>chr2\nAAC\n""")
    return genome_fa


@pytest.fixture
def bim_string():
    return """1 	rs11511647	0	1	A	T
1	rs3883674	0	4	C	G
1	rs12218882	0	6	T	T
2	rs12218882	0	2	T	T"""

@pytest.fixture
def bim_df(bim_string):

    return pd.read_table(StringIO(bim_string), sep=r"\s+", header=None,
                         names=["chromosome", "snp_name", "genetic_distance",
                                "position", "allele_1", "allele_2"],
                         dtype={"chromosome": str})
