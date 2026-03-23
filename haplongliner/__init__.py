"""
HapLongLINEr: A modular Python package for discovering and curating full-length (≥5 kb) young LINE-1 elements (L1HS, L1PA2, and potentially intact L1PA3) in haploid long-read human genome assemblies.

Modules:
    - rm_mode: RepeatMasker-based L1 discovery
    - mm_mode: minimap2-seeded, candidate-level RepeatMasker L1 discovery
    - sv_mode: SV-based (RepeatMasker-free) L1 discovery
    - sequence_retrieval_function: Pangenome-level L1 sequence retrieval

Authors:
    Lei Yang, Amanda Norseen, Rick McLaughlin
"""

from .rm_mode import *
from .mm_mode import *
from .sv_mode import *
from .sequence_retrieval_function import *
from .find_intact_orf import find_intact_orf
from .combine_table import combine_table
