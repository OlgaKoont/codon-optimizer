"""
Base specification classes for constraints and objectives.
"""

from abc import ABC, abstractmethod
from typing import Optional, List, Tuple
from dataclasses import dataclass


@dataclass
class Location:
    """Represents a location in a sequence."""
    start: int
    end: int
    strand: int = 1  # 1 for forward, -1 for reverse
    
    def __post_init__(self):
        """Validate location."""
        if self.start < 0 or self.end < 0:
            raise ValueError("Location start and end must be non-negative")
        if self.start >= self.end:
            raise ValueError("Location start must be less than end")
        if self.strand not in (1, -1):
            raise ValueError("Strand must be 1 or -1")
    
    @property
    def length(self) -> int:
        """Get length of location."""
        return self.end - self.start
    
    def overlaps(self, other: 'Location') -> bool:
        """Check if this location overlaps with another."""
        return not (self.end <= other.start or other.end <= self.start)
    
    def contains(self, other: 'Location') -> bool:
        """Check if this location contains another."""
        return self.start <= other.start and self.end >= other.end


class Specification(ABC):
    """
    Base class for all specifications (constraints and objectives).
    
    Inspired by DnaChisel's Specification system.
    """
    
    def __init__(self, 
                 location: Optional[Location] = None,
                 boost: float = 1.0,
                 priority: int = 0):
        """
        Initialize specification.
        
        Args:
            location: Location in sequence where specification applies (None = whole sequence)
            boost: Boost factor for objectives (relative importance)
            priority: Priority for solving (higher = solved first)
        """
        self.location = location
        self.boost = boost
        self.priority = priority
        self.best_possible_score: Optional[float] = None
        self.enforced_by_nucleotide_restrictions: bool = False
    
    @abstractmethod
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> 'SpecEvaluation':
        """
        Evaluate specification on a sequence.
        
        Args:
            sequence: DNA sequence
            protein_sequence: Protein sequence (if available)
            
        Returns:
            SpecEvaluation object
        """
        pass
    
    def localized(self, location: Location, problem=None, with_righthand: bool = True) -> Optional['Specification']:
        """
        Return a localized version of this specification for the given location.
        
        If the specification doesn't apply to this location, returns None.
        If it applies globally, returns self.
        Otherwise, returns a modified version scoped to the location.
        
        Args:
            location: Location to localize to
            problem: Optional problem instance for context
            with_righthand: If False, only consider left-hand side of location
            
        Returns:
            Localized specification or None
        """
        if self.location is None:
            # Global specification - return localized version
            return self.copy_with_changes(location=location)
        elif self.location.overlaps(location):
            # Overlaps - return intersection
            new_start = max(self.location.start, location.start)
            if with_righthand:
                new_end = min(self.location.end, location.end)
            else:
                new_end = min(self.location.end, location.start + (location.end - location.start) // 2)
            if new_start < new_end:
                return self.copy_with_changes(location=Location(new_start, new_end, self.location.strand))
        return None
    
    def copy_with_changes(self, **kwargs) -> 'Specification':
        """
        Return a copy of this specification with modified properties.
        
        Args:
            **kwargs: Properties to change
            
        Returns:
            New specification instance
        """
        import copy
        new_spec = copy.copy(self)
        for key, value in kwargs.items():
            setattr(new_spec, key, value)
        return new_spec
    
    def __repr__(self):
        loc_str = f"location={self.location}" if self.location else "location=global"
        return f"{self.__class__.__name__}({loc_str}, boost={self.boost}, priority={self.priority})"


class Constraint(Specification):
    """
    Base class for constraints (hard requirements that must be satisfied).
    """
    
    def __init__(self, **kwargs):
        """Initialize constraint."""
        super().__init__(**kwargs)
        self.is_focus: bool = False  # Focus constraint gets priority
    
    def passes(self, sequence: str, protein_sequence: Optional[str] = None) -> bool:
        """
        Check if constraint passes.
        
        Args:
            sequence: DNA sequence
            protein_sequence: Protein sequence (if available)
            
        Returns:
            True if constraint passes
        """
        evaluation = self.evaluate(sequence, protein_sequence)
        return evaluation.passes


class Objective(Specification):
    """
    Base class for objectives (soft goals to optimize).
    """
    
    def __init__(self, **kwargs):
        """Initialize objective."""
        super().__init__(**kwargs)
        self.optimize_passively: bool = False  # If True, only considered when optimizing other objectives


