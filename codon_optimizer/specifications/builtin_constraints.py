"""
Built-in constraint specifications.
"""

from typing import Optional, List, Tuple
from Bio.Seq import Seq
from .specification import Constraint, Location
from .evaluation import SpecEvaluation
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from codon_optimizer.analysis.gc_content import GCContentAnalyzer
from codon_optimizer.analysis.motifs import MotifDetector
from codon_optimizer.constraints.restriction_sites import RestrictionSiteManager
from codon_optimizer.constraints.protein_integrity import ProteinIntegrityChecker


class EnforceTranslation(Constraint):
    """
    Constraint to enforce that DNA sequence translates to a specific protein.
    """
    
    def __init__(self, protein_sequence: str, location: Optional[Location] = None, **kwargs):
        """
        Initialize translation constraint.
        
        Args:
            protein_sequence: Expected protein sequence
            location: Location where translation should be enforced
        """
        super().__init__(location=location, **kwargs)
        self.protein_sequence = protein_sequence
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate if sequence translates correctly."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
        else:
            seq_to_check = sequence
        
        # Ensure length is multiple of 3
        if len(seq_to_check) % 3 != 0:
            return SpecEvaluation(
                passes=False,
                score=0.0,
                message=f"Sequence length ({len(seq_to_check)}) is not multiple of 3",
                locations=[(self.location.start if self.location else 0, 
                           self.location.end if self.location else len(sequence))]
            )
        
        try:
            translated = str(Seq(seq_to_check).translate())
            if translated == self.protein_sequence:
                return SpecEvaluation(
                    passes=True,
                    score=1.0,
                    message="Translation matches expected protein sequence"
                )
            else:
                mismatches = sum(1 for a, b in zip(translated, self.protein_sequence) if a != b)
                return SpecEvaluation(
                    passes=False,
                    score=0.0,
                    message=f"Translation mismatch: {mismatches} differences",
                    locations=[(self.location.start if self.location else 0,
                               self.location.end if self.location else len(sequence))]
                )
        except Exception as e:
            return SpecEvaluation(
                passes=False,
                score=0.0,
                message=f"Translation error: {str(e)}"
            )


class EnforceGCContent(Constraint):
    """
    Constraint to enforce GC content within specified range.
    """
    
    def __init__(self, min_gc: float = 0.3, max_gc: float = 0.8, 
                 window: Optional[int] = None, location: Optional[Location] = None, **kwargs):
        """
        Initialize GC content constraint.
        
        Args:
            min_gc: Minimum GC content (0-1)
            max_gc: Maximum GC content (0-1)
            window: Window size for sliding window analysis (None = whole sequence)
            location: Location where constraint applies
        """
        super().__init__(location=location, **kwargs)
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.window = window
        self.gc_analyzer = GCContentAnalyzer(gc_min=min_gc, gc_max=max_gc, window_size=window or 50)
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate GC content constraint."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
            offset = self.location.start
        else:
            seq_to_check = sequence
            offset = 0
        
        violations = []
        overall_gc = (seq_to_check.count('G') + seq_to_check.count('C')) / len(seq_to_check) if seq_to_check else 0
        
        if self.window and len(seq_to_check) > self.window:
            # Check sliding windows
            for i in range(0, len(seq_to_check) - self.window + 1):
                window_seq = seq_to_check[i:i + self.window]
                window_gc = (window_seq.count('G') + window_seq.count('C')) / len(window_seq)
                if window_gc < self.min_gc or window_gc > self.max_gc:
                    violations.append((offset + i, offset + i + self.window))
        else:
            # Check overall
            if overall_gc < self.min_gc or overall_gc > self.max_gc:
                violations.append((offset, offset + len(seq_to_check)))
        
        if violations:
            score = max(0.0, 1.0 - len(violations) * 0.1)
            return SpecEvaluation(
                passes=False,
                score=score,
                message=f"GC content violation: {overall_gc:.3f} (target: {self.min_gc:.3f}-{self.max_gc:.3f})",
                locations=violations
            )
        else:
            return SpecEvaluation(
                passes=True,
                score=1.0,
                message=f"GC content OK: {overall_gc:.3f}"
            )


class AvoidPattern(Constraint):
    """
    Constraint to avoid specific DNA patterns (e.g., restriction sites).
    """
    
    def __init__(self, pattern: str, location: Optional[Location] = None, **kwargs):
        """
        Initialize pattern avoidance constraint.
        
        Args:
            pattern: DNA pattern to avoid (e.g., "BsaI_site", "GGATCC")
            location: Location where constraint applies
        """
        super().__init__(location=location, **kwargs)
        self.pattern = pattern.upper()
        # Map common restriction site names to sequences
        self.restriction_sites = {
            'BSAI_SITE': 'GGTCTC',
            'BSMBI_SITE': 'CGTCTC',
            'NOTI_SITE': 'GCGGCCGC',
            'ECORI_SITE': 'GAATTC',
            'BAMHI_SITE': 'GGATCC',
        }
        if pattern.upper() in self.restriction_sites:
            self.pattern = self.restriction_sites[pattern.upper()]
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate pattern avoidance."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
            offset = self.location.start
        else:
            seq_to_check = sequence
            offset = 0
        
        matches = []
        pattern_upper = self.pattern.upper()
        seq_upper = seq_to_check.upper()
        
        i = 0
        while i < len(seq_upper):
            pos = seq_upper.find(pattern_upper, i)
            if pos == -1:
                break
            matches.append((offset + pos, offset + pos + len(pattern_upper)))
            i = pos + 1
        
        if matches:
            score = max(0.0, 1.0 - len(matches) * 0.2)
            return SpecEvaluation(
                passes=False,
                score=score,
                message=f"Found {len(matches)} occurrences of pattern '{self.pattern}'",
                locations=matches
            )
        else:
            return SpecEvaluation(
                passes=True,
                score=1.0,
                message=f"No occurrences of pattern '{self.pattern}' found"
            )


class AvoidMotifs(Constraint):
    """
    Constraint to avoid problematic sequence motifs.
    """
    
    def __init__(self, motif_types: Optional[List[str]] = None, location: Optional[Location] = None, **kwargs):
        """
        Initialize motif avoidance constraint.
        
        Args:
            motif_types: List of motif types to avoid (None = all problematic motifs)
            location: Location where constraint applies
        """
        super().__init__(location=location, **kwargs)
        self.motif_types = motif_types or ['tata_box', 'cryptic_promoter', 'internal_rbs', 'polyA_signal']
        self.motif_detector = MotifDetector()
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate motif avoidance."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
            offset = self.location.start
        else:
            seq_to_check = sequence
            offset = 0
        
        all_motifs = self.motif_detector.find_motifs(seq_to_check)
        problematic_motifs = []
        
        for motif_type in self.motif_types:
            if motif_type in all_motifs:
                for start, end, match in all_motifs[motif_type]:
                    problematic_motifs.append((offset + start, offset + end))
        
        if problematic_motifs:
            # For long sequences, allow some motifs (soft constraint)
            # Score decreases with number of motifs, but doesn't require zero
            # For sequences > 400 bp, allow up to 5% of sequence length in motifs
            max_allowed = max(1, int(len(seq_to_check) * 0.05 / 10))  # Rough estimate
            score = max(0.0, 1.0 - (len(problematic_motifs) - max_allowed) * 0.1) if len(problematic_motifs) > max_allowed else 0.8
            passes = len(problematic_motifs) <= max_allowed
            
            return SpecEvaluation(
                passes=passes,
                score=max(0.1, score),  # Never return 0.0 to allow gradient search
                message=f"Found {len(problematic_motifs)} problematic motifs (max allowed: {max_allowed})",
                locations=problematic_motifs
            )
        else:
            return SpecEvaluation(
                passes=True,
                score=1.0,
                message="No problematic motifs found"
            )


class EnforceRestrictionSites(Constraint):
    """
    Constraint to enforce absence/presence of restriction sites.
    """
    
    def __init__(self, sites_to_remove: List[str], sites_to_keep: Optional[List[str]] = None,
                 location: Optional[Location] = None, **kwargs):
        """
        Initialize restriction site constraint.
        
        Args:
            sites_to_remove: List of restriction sites that must be absent
            sites_to_keep: List of restriction sites that must be present (optional)
            location: Location where constraint applies
        """
        super().__init__(location=location, **kwargs)
        self.sites_to_remove = sites_to_remove
        self.sites_to_keep = sites_to_keep or []
        self.restriction_manager = RestrictionSiteManager(
            sites_to_remove=sites_to_remove,
            sites_to_keep=sites_to_keep
        )
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate restriction site constraint."""
        violations = self.restriction_manager.check_violations(sequence)
        
        if violations['has_violations']:
            locations = [(v['start'], v['end']) for v in violations.get('violations', [])]
            return SpecEvaluation(
                passes=False,
                score=0.5,
                message=f"Restriction site violations: {len(locations)}",
                locations=locations
            )
        else:
            return SpecEvaluation(
                passes=True,
                score=1.0,
                message="All restriction site constraints satisfied"
            )


class EnforceProteinIntegrity(Constraint):
    """
    Constraint to enforce that optimized sequence encodes the same protein.
    """
    
    def __init__(self, original_protein: str, location: Optional[Location] = None, **kwargs):
        """
        Initialize protein integrity constraint.
        
        Args:
            original_protein: Original protein sequence that must be preserved
            location: Location where constraint applies
        """
        super().__init__(location=location, **kwargs)
        self.original_protein = original_protein
        self.integrity_checker = ProteinIntegrityChecker()
    
    def evaluate(self, sequence: str, protein_sequence: Optional[str] = None) -> SpecEvaluation:
        """Evaluate protein integrity."""
        if self.location:
            seq_to_check = sequence[self.location.start:self.location.end]
        else:
            seq_to_check = sequence
        
        try:
            optimized_protein = str(Seq(seq_to_check).translate())
            violations = self.integrity_checker.check_integrity_violations(
                self.original_protein, optimized_protein
            )
            
            if violations['is_valid']:
                return SpecEvaluation(
                    passes=True,
                    score=1.0,
                    message="Protein integrity preserved"
                )
            else:
                return SpecEvaluation(
                    passes=False,
                    score=0.0,
                    message=f"Protein integrity violation: {violations.get('errors', [])}",
                    details=violations
                )
        except Exception as e:
            return SpecEvaluation(
                passes=False,
                score=0.0,
                message=f"Translation error: {str(e)}"
            )

