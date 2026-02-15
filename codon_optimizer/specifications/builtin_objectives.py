"""
Built-in objective specifications.
"""

from typing import Optional
from .specification import Objective, Location
from .evaluation import SpecEvaluation
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from codon_optimizer.analysis.gc_content import GCContentAnalyzer
from codon_optimizer.analysis.motifs import MotifDetector
from codon_optimizer.analysis.codon_pairs import CodonPairAnalyzer
from codon_optimizer.analysis.mrna_structure import MRNAStructureAnalyzer


class MaximizeCAI(Objective):
    """
    Objective to maximize Codon Adaptation Index (CAI).
    """
    
    def __init__(self, host: str = "CHO", target_min: float = 0.7, target_max: float = 0.9,
                 location: Optional[Location] = None, **kwargs):
        """
        Initialize CAI maximization objective.
        
        Args:
            host: Host organism
            target_min: Minimum target CAI
            target_max: Maximum target CAI (to avoid over-optimization)
            location: Location where objective applies
        """
        super().__init__(location=location, **kwargs)
        self.host = host
        self.target_min = target_min
        self.target_max = target_max
        self.codon_analyzer = CodonUsageAnalyzer(host=host)
        self.best_possible_score = 1.0
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate CAI score."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
        else:
            seq_to_check = sequence
        
        cai = self.codon_analyzer.calculate_cai(seq_to_check)
        
        # Score based on how close to target range
        if self.target_min <= cai <= self.target_max:
            # In target range - normalize to 0-1
            score = (cai - self.target_min) / (self.target_max - self.target_min)
        elif cai < self.target_min:
            # Below target - penalty
            score = cai / self.target_min
        else:
            # Above target - slight penalty for over-optimization
            score = 1.0 - (cai - self.target_max) / (1.0 - self.target_max)
            score = max(0.0, score)
        
        return SpecEvaluation(
            passes=score > 0.5,
            score=score,
            message=f"CAI: {cai:.3f} (target: {self.target_min:.3f}-{self.target_max:.3f})"
        )


class OptimizeCodonUsage(Objective):
    """
    Objective to optimize codon usage (avoid tRNA depletion).
    """
    
    def __init__(self, host: str = "CHO", max_codon_usage: float = 0.5,
                 location: Optional[Location] = None, **kwargs):
        """
        Initialize codon usage optimization objective.
        
        Args:
            host: Host organism
            max_codon_usage: Maximum allowed usage of single codon
            location: Location where objective applies
        """
        super().__init__(location=location, **kwargs)
        self.host = host
        self.max_codon_usage = max_codon_usage
        self.codon_analyzer = CodonUsageAnalyzer(host=host, max_codon_usage=max_codon_usage)
        self.best_possible_score = 1.0
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate codon usage optimization."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
        else:
            seq_to_check = sequence
        
        has_depletion, depletion_info = self.codon_analyzer.check_tRNA_depletion(seq_to_check)
        
        if has_depletion:
            # Penalty for tRNA depletion risk
            max_usage = max(depletion_info.get('usage_frequencies', {}).values(), default=0.0)
            score = max(0.0, 1.0 - (max_usage - self.max_codon_usage) / (1.0 - self.max_codon_usage))
            return SpecEvaluation(
                passes=False,
                score=score,
                message=f"tRNA depletion risk: max codon usage = {max_usage:.3f}"
            )
        else:
            return SpecEvaluation(
                passes=True,
                score=1.0,
                message="Codon usage balanced (no tRNA depletion risk)"
            )


class OptimizeGCContent(Objective):
    """
    Objective to optimize GC content (keep within target range).
    """
    
    def __init__(self, target_gc: float = 0.5, tolerance: float = 0.1,
                 location: Optional[Location] = None, **kwargs):
        """
        Initialize GC content optimization objective.
        
        Args:
            target_gc: Target GC content (0-1)
            tolerance: Acceptable deviation from target
            location: Location where objective applies
        """
        super().__init__(location=location, **kwargs)
        self.target_gc = target_gc
        self.tolerance = tolerance
        self.gc_analyzer = GCContentAnalyzer(
            gc_min=target_gc - tolerance,
            gc_max=target_gc + tolerance,
            window_size=50
        )
        self.best_possible_score = 1.0
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate GC content optimization."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
        else:
            seq_to_check = sequence
        
        score = self.gc_analyzer.calculate_gc_score(seq_to_check)
        overall_gc = (seq_to_check.count('G') + seq_to_check.count('C')) / len(seq_to_check) if seq_to_check else 0
        
        return SpecEvaluation(
            passes=score > 0.7,
            score=score,
            message=f"GC content: {overall_gc:.3f} (target: {self.target_gc:.3f} ± {self.tolerance:.3f})"
        )


class MinimizeMotifs(Objective):
    """
    Objective to minimize problematic sequence motifs.
    """
    
    def __init__(self, motif_types: Optional[list] = None, location: Optional[Location] = None, **kwargs):
        """
        Initialize motif minimization objective.
        
        Args:
            motif_types: List of motif types to minimize (None = all)
            location: Location where objective applies
        """
        super().__init__(location=location, **kwargs)
        self.motif_types = motif_types
        self.motif_detector = MotifDetector()
        self.best_possible_score = 1.0
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate motif minimization."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
        else:
            seq_to_check = sequence
        
        all_motifs = self.motif_detector.find_motifs(seq_to_check)
        total_motifs = sum(len(motifs) for motifs in all_motifs.values())
        
        # Score decreases with number of motifs
        score = max(0.0, 1.0 / (1.0 + total_motifs * 0.1))
        
        return SpecEvaluation(
            passes=total_motifs == 0,
            score=score,
            message=f"Found {total_motifs} problematic motifs"
        )


class OptimizeCodonPairs(Objective):
    """
    Objective to optimize codon pair usage.
    """
    
    def __init__(self, host: str = "CHO", location: Optional[Location] = None, **kwargs):
        """
        Initialize codon pair optimization objective.
        
        Args:
            host: Host organism
            location: Location where objective applies
        """
        super().__init__(location=location, **kwargs)
        self.host = host
        self.pair_analyzer = CodonPairAnalyzer(host=host)
        self.best_possible_score = 1.0
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate codon pair optimization."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
        else:
            seq_to_check = sequence
        
        score = self.pair_analyzer.calculate_pair_score(seq_to_check)
        
        return SpecEvaluation(
            passes=score > 0.7,
            score=score,
            message=f"Codon pair optimization score: {score:.3f}"
        )


class OptimizeMRNAStructure(Objective):
    """
    Objective to optimize mRNA secondary structure (minimize stable structures in 5' region).
    """
    
    def __init__(self, max_5prime_dg: float = -10.0, location: Optional[Location] = None, **kwargs):
        """
        Initialize mRNA structure optimization objective.
        
        Args:
            max_5prime_dg: Maximum ΔG for 5' region (kcal/mol)
            location: Location where objective applies
        """
        super().__init__(location=location, **kwargs)
        self.max_5prime_dg = max_5prime_dg
        self.structure_analyzer = MRNAStructureAnalyzer(
            max_5prime_dg=max_5prime_dg,
            structure_region_start=-50,
            structure_region_end=100
        )
        self.best_possible_score = 1.0
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate mRNA structure optimization."""
        start_codon_pos = sequence.find('ATG') if 'ATG' in sequence else 0
        score = self.structure_analyzer.calculate_structure_score(sequence, start_codon_pos)
        
        return SpecEvaluation(
            passes=score > 0.7,
            score=score,
            message=f"mRNA structure score: {score:.3f}"
        )

