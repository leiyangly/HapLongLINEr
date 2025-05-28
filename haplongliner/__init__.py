"""
HapLongLINEr: A modular Python package for discovering and curating full-length (≥5 kb) young LINE-1 elements (L1HS, L1PA2, and potentially intact L1PA3) in haploid long-read human genome assemblies.

Modules:
    - module1_RM: RepeatMasker-based L1 discovery
    - module2_SV: SV-based (RepeatMasker-free) L1 discovery
    - module3_DB: Pangenome-level L1 sequence repository

Authors:
    Lei Yang, Amanda Norseen, Rick McLaughlin
"""

from .module1_RM import *
from .module2_SV import *
from .module3_DB import *