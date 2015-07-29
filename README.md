#snpflip

snpflip outputs each reverse strand and ambiguous strand SNP in your GWAS data.

##Main use cases

- **Generate a list of SNPs not on the reference strand:** Many bioinformatics applications require the SNPs to be on the reference strand. By running snpflip on your GWAS data you'll find the SNPs that need to be flipped.

- **Quality controlling your GWAS data and pruning it:** If your GWAS data contains ambiguous SNPs, these might ruin the usefulness of your data for certain purposes (e.g., imputation). By running snpflip you can get the names of the ambiguous SNPs and remove them with Plink.

##Install
`pip install snpflip`

##Usage

```
snpflip

Report reverse strand and ambiguous strand SNPs.
(Visit github.com/endrebak/snp-flip for examples and help.)

Usage:
    snpflip --fasta-genome=FA --bim-file=BIM [--output-prefix=PRE]
    snpflip --help

Arguments:
    -f FA --fasta-genome=FA    fasta genome
    -b BIM --bim-file=BIM      plink bed file (extended variant information file)

Options:
    -h --help                   show this message
    -o PRE --output-prefix=PRE  the prefix of the output-files
    (./snpflip_output/<bim_basename> by default)
```

##Output

The output files are:
- `<prefix>.reverse` - The SNPs that are on the reverse strand.
- `<prefix>.ambiguous` - The SNPs that cannot be assigned to a strand.
- `<prefix>.annotated_bim` - Strand annotated bim table.


which means the snpflip output files can be used as input to Plink. This is convenient if you want to remove the bad or ambiguous SNPs or to flip the SNPs that are on the reverse strand.

snpflip is not able to find out whether an insertion/deletion or structural variant is on the correct strand. If you find reverse strand SNPs in your dataset, this is an indication that you cannot trust the strandedness of your data. Hence, you might want to remove all structural variants and insertions/deletions. This can be done by giving plink the flag `--snps-only no-DI`.

####Example usage with plink:

```python2 snipflip.py -b snp_data.bim -rgf genome.fa -op snpflip_output```

```plink --file snp_data --flip snpflip_output.reverse --make-bed```


##Extended example

This example uses the example files in the example_files catalog.

```
endrebak@tang:~/biocore/snpflip$ cat example_files/reference_genome.fa
>chr1
ACT
>chr2
CCC
endrebak@tang:~/biocore/snpflip$ cat example_files/example.bim
1 snp1 0 1 A C
1 snp2 0 2 A T
1 snp3 0 3 A G
2 snp4 0 1 A G
2 esv5 0 2 AA G
endrebak@tang:~/biocore/snpflip$ python2 snpflip.py -b example_files/example.bim -rgf \
 example_files/reference_genome.fa -op extended_example
Loading reference genome file example_files/reference_genome.fa
***
Results for example_files/example.bim
***
                          Reverse strand      Forward strand       Ambiguous snp
         Monomorphic                   0                   0                   0
         Polymorphic                   2                   1                   1
  Structural variant                  NA                  NA                   1
         Invalid snp                  NA                  NA                   0
Output files written with prefix extended_example
endrebak@tang:~/biocore/snpflip$ cat extended_example.ambiguous
snp2    1    2    A    T    C
endrebak@tang:~/biocore/snpflip$ cat extended_example.reverse
snp3    1    3    A    G    T
snp4	2    1    A    G    C
```

##Requirements

The only third-party module requirement is [Biopython](http://biopython.org/). It is only used in a few lines (~5) to do reference genome lookups and can easily be changed to use the library of your choice.

To allow snpflip to map the chromosomes in the .bim and .fa files, the chromosome names need to be 1, 2, .., X, Y in the bim file and chr1, chr2, ..., chrX, chrY in the fasta file. (This is the standard convention.) See the files in the example_files catalog for examples.

##FAQ

####My SNP-data is not in Plink .bed+.bim+.fam format. What do I do?

Use [Plink](https://www.cog-genomics.org/plink2/data) to convert your files to the binary plink format. Example:

`plink --vcf yourfile.vcf --make-bed --out your_prefix`

####How do you decide the SNP type and what strand it came from?

Since the reference genome is based on the forward strand, finding an 'A' in the reference genome and the bim-file for said SNP tells you that it is from the forward strand. If the SNP were a 'T' it would have been a reverse strand SNP.

####Why do I get a 'KeyError'?

You need to make sure that your bim and fasta files follow the chromosome naming conventions outlined in the "Requirements" section above.


##Issues
See https://github.com/endrebak/snipflip/issues

No issue is too small. I will entertain all feature requests. Please ask support questions by raising an issue.

##Todo
Add unit-tests.
