
import pytest

import pandas as pd

from ebs.imports import StringIO

from snp_flip.bim_file import read_bim_file

def describe_get_snp_data():

    def returns_dataframe(bim_file, expected_df):

        actual_df = read_bim_file(str(bim_file))
        print(actual_df)
        print(expected_df)
        assert actual_df.equals(expected_df)


@pytest.fixture
def bim_file(tmpdir, bim_string):
    f = tmpdir.join("test")
    f.write(bim_string)
    return f


@pytest.fixture
def bim_string():
    return """21	rs11511647	0	26765	A	T
X	rs3883674	0	32380	C	G
X	rs12218882	0	48172	T	T
9	rs10904045	0	48426	A	T
9	rs10751931	0	49949	C	T
8	rs11252127	0	52087	A	C
10	rs12775203	0	52277	A	A
8	rs12255619	0	52481	G	T"""

@pytest.fixture
def expected_df(bim_string):

    return pd.read_table(StringIO(bim_string), sep=r"\s+", header=None,
                         names=["chromosome", "snp_name", "genetic_distance",
                                "position", "allele_1", "allele_2"])
