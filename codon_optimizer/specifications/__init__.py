"""
Specification system for codon optimization.
Inspired by DnaChisel architecture.
"""

from .specification import Specification, Constraint, Objective
from .evaluation import SpecEvaluation, SpecEvaluations
from .builtin_constraints import (
    EnforceTranslation,
    EnforceGCContent,
    AvoidPattern,
    AvoidMotifs,
    EnforceRestrictionSites,
    EnforceProteinIntegrity,
)
from .builtin_objectives import (
    MaximizeCAI,
    OptimizeCodonUsage,
    OptimizeGCContent,
    MinimizeMotifs,
    OptimizeCodonPairs,
    OptimizeMRNAStructure,
)

__all__ = [
    'Specification',
    'Constraint',
    'Objective',
    'SpecEvaluation',
    'SpecEvaluations',
    'EnforceTranslation',
    'EnforceGCContent',
    'AvoidPattern',
    'AvoidMotifs',
    'EnforceRestrictionSites',
    'EnforceProteinIntegrity',
    'MaximizeCAI',
    'OptimizeCodonUsage',
    'OptimizeGCContent',
    'MinimizeMotifs',
    'OptimizeCodonPairs',
    'OptimizeMRNAStructure',
]

