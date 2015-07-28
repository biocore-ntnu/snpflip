
import pytest

import pandas as pd


from ebs.imports import StringIO

from snp_flip.table import create_snp_table

def test_create_snp_table(bim_file, fa_file, expected_result):

    actual_result = create_snp_table(bim_file, fa_file)
    print(actual_result)
    print(expected_result)

    assert actual_result.equals(expected_result)

@pytest.fixture
def fa_file(tmpdir):

    genome_fa = tmpdir.join("genome.fa")
    genome_fa.write(u""">chr1\nACGTTCCG\n>chr2\nAAC\n""")
    return str(genome_fa)


@pytest.fixture
def bim_file(bim_string, tmpdir):

    bim = tmpdir.join("b.bim")
    bim.write(bim_string)
    return str(bim)


@pytest.fixture
def bim_string():
    return """1	rs11252127	0	0	A	C
1	rs12255619	0	3	G	T
2	rs10751931	0	2	C	T"""


@pytest.fixture
def expected_result():
    return pd.read_table(StringIO("""1	rs11252127	0	0	A	C   A   T
1	rs12255619	0	3	G	T   T   A
2	rs10751931	0	2	C	T   C   G"""),
                        sep=r"\s+", header=None,
                        names=["chromosome", "snp_name", "genetic_distance",
                               "position", "allele_1", "allele_2", "reference",
                               "reference_rev"], dtype={"chromosome": str})

# def compute_strand(args):
