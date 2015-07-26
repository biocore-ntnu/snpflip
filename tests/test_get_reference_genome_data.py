import pytest



from snp_flip.reference_genome import get_reference_genome_data

def describe_get_reference_genome_data():

    def test_get_data(snp_positions, genome_fasta):
        snp_nucleotides = get_reference_genome_data(snp_positions, "chr1", str(genome_fasta))
        assert snp_nucleotides == [u"A", u"T", u"C"]



@pytest.fixture
def genome_fasta(tmpdir):

    genome_fa = tmpdir.join("genome.fa")
    genome_fa.write(u""">chr1\nACGTTCCG\n>chr2\nAAC\n""")
    return genome_fa

@pytest.fixture
def snp_positions():
    return [0, 3, 5]
