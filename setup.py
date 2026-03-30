from setuptools import setup, find_packages

setup(
    name="haplongliner",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "haplongliner": [
            "data/L1rp.fa",
            "data/L1rpORF12p.fa",
            "data/L1rpORF12p.fa.pdb",
            "data/L1rpORF12p.fa.phr",
            "data/L1rpORF12p.fa.pin",
            "data/L1rpORF12p.fa.pot",
            "data/L1rpORF12p.fa.psq",
            "data/L1rpORF12p.fa.ptf",
            "data/L1rpORF12p.fa.pto",
            "data/HPRC_L1_hs1_v2_v2fl_status.bed",
            "data/HPRC_L1_seq_by_site_v2.zip",
            "data/HPRC_L1_seq_by_site_v2fl.zip",
        ]
    },
    python_requires=">=3.10",
    install_requires=[
        "biopython",
        "pyfaidx",
        "edlib",
        "tqdm",
    ],
    entry_points={
        "console_scripts": [
            "haplongliner=haplongliner.cli:main",
            "liftover_paf=haplongliner.liftover_paf:main",
        ]
    },
)
