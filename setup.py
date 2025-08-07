from setuptools import setup, find_packages

setup(
    name="haplongliner",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
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
