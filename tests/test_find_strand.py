

import pytest

import pandas as pd

from snp_flip.table import find_strand

from ebs.imports import StringIO


def test_find_strand(input_df, expected_result):

    actual_df = find_strand(input_df)
    print(actual_df)
    print(expected_result)
    assert actual_df.equals(expected_result)



@pytest.fixture
def input_df():
    return pd.read_table(StringIO("""
1	rs11252127	0	0	A	C   A   T
1	rs12255619	0	3	G	T   T   A
2	rs10751931	0	2	C	T   C   G
3	rs11252128	0	0	G	C   A   T
3	rs12255620	0	3	A	T   T   A
3	rs10751932	0	2	G	T   C   G"""),
                        sep=r"\s+", header=None,
                        names=["chromosome", "snp_name", "genetic_distance",
                               "position", "allele_1", "allele_2", "reference",
                               "reference_rev"],
                         dtype={"chromosome": str})

@pytest.fixture
def expected_result():
    return pd.read_table(StringIO("""
1	rs11252127	0	0	A	C   A   T   forward
1	rs12255619	0	3	G	T   T   A   forward
2	rs10751931	0	2	C	T   C   G   forward
3	rs11252128	0	0	G	C   A   T   ambiguous
3	rs12255620	0	3	A	T   T   A   ambiguous
3	rs10751932	0	2	G	T   C   G   reverse"""),
                        sep=r"\s+", header=None,
                        names=["chromosome", "snp_name", "genetic_distance",
                               "position", "allele_1", "allele_2", "reference",
                               "reference_rev", "strand"],
                         dtype={"chromosome": str})
