"""
HapLongLINEr: A modular Python package for discovering and curating full-length (≥5 kb) young LINE-1 elements (L1HS, L1PA2, and potentially intact L1PA3) in haploid long-read human genome assemblies.

Modules:
    - rm_module: RepeatMasker-based L1 discovery
    - sv_module: SV-based (RepeatMasker-free) L1 discovery
    - module3_DB: Pangenome-level L1 sequence repository

Authors:
    Lei Yang, Amanda Norseen, Rick McLaughlin
"""

from .rm_module import *
from .sv_module import *
from .module3_DB import *
from .find_intact_orf import find_intact_orf
from .combine_table import combine_table
