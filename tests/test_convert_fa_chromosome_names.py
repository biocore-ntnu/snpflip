
import pytest

from snp_flip.reference_genome import convert_fa_chromosome_names

def test_convert_fa_chromosome_names(ucsc_dict, expected_result):

    fasta_dict_converted = convert_fa_chromosome_names(ucsc_dict)
    assert fasta_dict_converted == expected_result


@pytest.fixture
def ucsc_dict():
    return {"chr22": "", "chrX": "", "chrY": "", "chrM": "",
            "chr19_GL949747v2_alt": ""}


@pytest.fixture
def onekg_dict():
    return {"22": "", "X": "", "Y": "", "MT": "",
            "GL000207.1 dna:supercontig supercontig::GL000207.1:1:4262:1": ""}


@pytest.fixture
def expected_result():
    return {"22": "", "X": "", "Y": "", "M": ""}

# >^(?:chr)?(\d+|X|Y|M)(?:T)?\s+22 dna:chromosome chromosome:GRCh37:22:1:51304566:1
# >X dna:chromosome chromosome:GRCh37:X:1:155270560:1
# >Y dna:chromosome chromosome:GRCh37:Y:2649521:59034049:1
# >MT gi|251831106|ref|NC_012920.1| Homo sapiens mitochondrion, complete genome
# >GL000207.1 dna:supercontig supercontig::GL000207.1:1:4262:1
# >GL000226.1 dna:supercontig supercontig::GL000226.1:1:15008:1
# >GL000229.1 dna:supercontig supercontig::GL000229.1:1:19913:1
# >GL000231.1 dna:supercontig supercontig::GL000231.1:1:27386:1
# >GL000210.1 dna:supercontig supercontig::GL000210.1:1:27682:1
# >GL000239.1 dna:supercontig supercontig::GL000239.1:1:33824:1
# >GL000235.1 dna:supercontig supercontig::GL000235.1:1:34474:1
# >GL000201.1 dna:supercontig supercontig::GL000201.1:1:36148:1

# # >chrX
# # >chrY
# # >chrM
# >chr19_GL949747v2_alt
# >chr22_KB663609v1_alt
# >chrX_KI270913v1_alt
