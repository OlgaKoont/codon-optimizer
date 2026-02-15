"""
Evaluation results for specifications.
"""

from dataclasses import dataclass
from typing import Optional, List, Tuple, Dict
from .specification import Specification


@dataclass
class SpecEvaluation:
    """
    Result of evaluating a specification on a sequence.
    
    Attributes:
        passes: Whether constraint passes (for constraints) or score > 0 (for objectives)
        score: Score of the specification (0-1, higher is better)
        message: Human-readable message about the evaluation
        locations: List of (start, end) tuples where issues occur (for constraints)
        details: Additional details dictionary
    """
    passes: bool
    score: float
    message: str = ""
    locations: List[Tuple[int, int]] = None
    details: Dict = None
    
    def __post_init__(self):
        """Initialize default values."""
        if self.locations is None:
            self.locations = []
        if self.details is None:
            self.details = {}
    
    def __repr__(self):
        status = "PASS" if self.passes else "FAIL"
        return f"SpecEvaluation({status}, score={self.score:.3f}, {self.message})"


class SpecEvaluations:
    """Collection of specification evaluations."""
    
    def __init__(self, evaluations: List[Tuple[Specification, SpecEvaluation]]):
        """
        Initialize evaluations collection.
        
        Args:
            evaluations: List of (specification, evaluation) tuples
        """
        self.evaluations = evaluations
    
    def filter(self, status: str) -> 'SpecEvaluations':
        """
        Filter evaluations by status.
        
        Args:
            status: 'passing', 'failing', or 'all'
            
        Returns:
            Filtered evaluations
        """
        if status == 'passing':
            filtered = [(spec, eval_) for spec, eval_ in self.evaluations if eval_.passes]
        elif status == 'failing':
            filtered = [(spec, eval_) for spec, eval_ in self.evaluations if not eval_.passes]
        else:
            filtered = self.evaluations
        return SpecEvaluations(filtered)
    
    def scores_sum(self) -> float:
        """Sum of all scores (weighted by boost)."""
        return sum(spec.boost * eval_.score for spec, eval_ in self.evaluations)
    
    def to_text(self) -> str:
        """Convert to human-readable text."""
        lines = []
        for spec, eval_ in self.evaluations:
            status = "✓" if eval_.passes else "✗"
            lines.append(f"{status} {spec.__class__.__name__}: {eval_.message} (score: {eval_.score:.3f})")
        return "\n".join(lines)


