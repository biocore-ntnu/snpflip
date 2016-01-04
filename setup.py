from setuptools import setup

from snp_flip.version import __version__

setup(
    name = "snpflip",
    packages = ["snp_flip"],
    scripts = ["bin/snpflip"],
    version = __version__,
    description = "Report reverse and ambiguous strand SNPs.",
    author = "Endre Bakken Stovner",
    author_email = "endrebak@stud.ntnu.no",
    url = "http://github.com/endrebak/snp-flip",
    keywords = ["snp", "gwas", "strand"],
    license = ["GPL-3.0"],
    install_requires = ["pandas>=0.16", "ebs", "docopt", "pyfaidx", "natsort"],
    classifiers = [
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"],
    long_description = ("Report reverse and ambiguous strand SNPs.\n"
                        "See the URL for examples and docs.")
)
