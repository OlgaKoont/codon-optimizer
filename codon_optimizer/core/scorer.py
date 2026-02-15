"""
Multi-criteria scoring function for codon optimization.
"""

from typing import Dict, Optional
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from codon_optimizer.analysis.gc_content import GCContentAnalyzer
from codon_optimizer.analysis.codon_pairs import CodonPairAnalyzer
from codon_optimizer.analysis.mrna_structure import MRNAStructureAnalyzer
from codon_optimizer.analysis.motifs import MotifDetector
from codon_optimizer.constraints.restriction_sites import RestrictionSiteManager
from codon_optimizer.constraints.protein_integrity import ProteinIntegrityChecker
from codon_optimizer.constraints.expression_cassette import ExpressionCassetteChecker
from Bio.Seq import Seq


class MultiCriteriaScorer:
    """Multi-criteria scorer for codon optimization."""
    
    def __init__(self, config: OptimizationConfig):
        """
        Initialize scorer with configuration.
        
        Args:
            config: Optimization configuration
        """
        self.config = config
        
        # Initialize analyzers
        self.codon_analyzer = CodonUsageAnalyzer(
            host=config.host,
            max_codon_usage=config.max_codon_usage
        )
        self.gc_analyzer = GCContentAnalyzer(
            gc_min=config.gc_min,
            gc_max=config.gc_max,
            window_size=config.gc_window_size
        )
        self.pair_analyzer = CodonPairAnalyzer(host=config.host)
        self.structure_analyzer = MRNAStructureAnalyzer(
            max_5prime_dg=config.max_5prime_dg,
            structure_region_start=config.structure_region_start,
            structure_region_end=config.structure_region_end
        )
        self.motif_detector = MotifDetector()
        self.restriction_manager = RestrictionSiteManager(
            sites_to_remove=config.restriction_sites_to_remove,
            sites_to_keep=config.restriction_sites_to_keep
        )
        self.integrity_checker = ProteinIntegrityChecker()
        self.cassette_checker = ExpressionCassetteChecker(
            kozak_sequence=config.kozak_sequence
        )
    
    def calculate_cai_score(self, dna_sequence: str) -> float:
        """
        Calculate CAI-based score.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Score (0-1)
        """
        cai = self.codon_analyzer.calculate_cai(dna_sequence)
        
        # Check if CAI is in target range
        if self.config.target_cai_min <= cai <= self.config.target_cai_max:
            # Normalize to 0-1 within target range
            score = (cai - self.config.target_cai_min) / (
                self.config.target_cai_max - self.config.target_cai_min
            )
        elif cai < self.config.target_cai_min:
            # Below target - penalty
            score = cai / self.config.target_cai_min
        else:
            # Above target - slight penalty for over-optimization
            score = 1.0 - (cai - self.config.target_cai_max) / (1.0 - self.config.target_cai_max)
            score = max(0.0, score)
        
        # Check for tRNA depletion
        has_depletion, _ = self.codon_analyzer.check_tRNA_depletion(dna_sequence)
        if has_depletion:
            score *= 0.7  # Penalty for depletion risk
        
        return max(0.0, min(1.0, score))
    
    def calculate_gc_score(self, dna_sequence: str) -> float:
        """
        Calculate GC content score.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Score (0-1)
        """
        return self.gc_analyzer.calculate_gc_score(dna_sequence)
    
    def calculate_pair_score(self, dna_sequence: str) -> float:
        """
        Calculate codon pair preference score.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Score (0-1)
        """
        return self.pair_analyzer.calculate_pair_score(dna_sequence)
    
    def calculate_structure_score(self, dna_sequence: str, start_codon_pos: int = 0) -> float:
        """
        Calculate mRNA structure score.
        
        Args:
            dna_sequence: DNA sequence
            start_codon_pos: Position of start codon
            
        Returns:
            Score (0-1)
        """
        return self.structure_analyzer.calculate_structure_score(dna_sequence, start_codon_pos)
    
    def calculate_motif_score(self, dna_sequence: str) -> float:
        """
        Calculate motif penalty score (inverted - higher is better).
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Score (0-1)
        """
        penalty = self.motif_detector.calculate_motif_penalty(dna_sequence)
        return 1.0 - penalty  # Invert penalty to score
    
    def calculate_cloning_score(self, dna_sequence: str) -> float:
        """
        Calculate cloning compatibility score.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Score (0-1)
        """
        return self.restriction_manager.calculate_cloning_score(dna_sequence)
    
    def calculate_total_score(self,
                             dna_sequence: str,
                             original_protein: Optional[str] = None,
                             start_codon_pos: int = 0) -> Dict[str, float]:
        """
        Calculate total multi-criteria score.
        
        Args:
            dna_sequence: DNA sequence
            original_protein: Original protein sequence (for validation)
            start_codon_pos: Position of start codon
            
        Returns:
            Dictionary with individual and total scores
        """
        # Validate protein integrity if original provided
        if original_protein:
            optimized_protein = str(Seq(dna_sequence).translate())
            violations = self.integrity_checker.check_integrity_violations(
                original_protein, optimized_protein
            )
            if not violations['is_valid']:
                # Severe penalty for integrity violations
                return {
                    'cai_score': 0.0,
                    'gc_score': 0.0,
                    'pair_score': 0.0,
                    'structure_score': 0.0,
                    'motif_score': 0.0,
                    'cloning_score': 0.0,
                    'total_score': 0.0,
                    'integrity_violation': True
                }
        
        # Calculate individual scores
        cai_score = self.calculate_cai_score(dna_sequence)
        gc_score = self.calculate_gc_score(dna_sequence)
        pair_score = self.calculate_pair_score(dna_sequence)
        structure_score = self.calculate_structure_score(dna_sequence, start_codon_pos)
        motif_score = self.calculate_motif_score(dna_sequence)
        cloning_score = self.calculate_cloning_score(dna_sequence)
        
        # Calculate weighted total
        total_score = (
            self.config.weight_cai * cai_score +
            self.config.weight_gc * gc_score +
            self.config.weight_codon_pairs * pair_score +
            self.config.weight_mrna_structure * structure_score +
            self.config.weight_motifs * motif_score +
            self.config.weight_cloning * cloning_score
        )
        
        return {
            'cai_score': cai_score,
            'gc_score': gc_score,
            'pair_score': pair_score,
            'structure_score': structure_score,
            'motif_score': motif_score,
            'cloning_score': cloning_score,
            'total_score': total_score,
            'integrity_violation': False
        }
    
    def get_detailed_analysis(self,
                              dna_sequence: str,
                              original_protein: Optional[str] = None) -> Dict:
        """
        Get detailed analysis of sequence.
        
        Args:
            dna_sequence: DNA sequence
            original_protein: Original protein sequence
            
        Returns:
            Dictionary with detailed analysis
        """
        start_codon_pos = self.cassette_checker.find_start_codon(dna_sequence) or 0
        
        scores = self.calculate_total_score(dna_sequence, original_protein, start_codon_pos)
        
        # Get detailed statistics
        cai = self.codon_analyzer.calculate_cai(dna_sequence)
        gc_stats = self.gc_analyzer.get_gc_statistics(dna_sequence)
        pair_stats = self.pair_analyzer.get_pair_statistics(dna_sequence)
        structure_stats = self.structure_analyzer.get_structure_statistics(
            dna_sequence, start_codon_pos
        )
        motifs = self.motif_detector.get_all_problematic_motifs(dna_sequence)
        restriction_violations = self.restriction_manager.check_violations(dna_sequence)
        cassette_stats = self.cassette_checker.get_cassette_statistics(dna_sequence)
        
        analysis = {
            'scores': scores,
            'cai': cai,
            'gc_statistics': gc_stats,
            'pair_statistics': pair_stats,
            'structure_statistics': structure_stats,
            'motifs': motifs,
            'restriction_sites': restriction_violations,
            'expression_cassette': cassette_stats,
        }
        
        if original_protein:
            optimized_protein = str(Seq(dna_sequence).translate())
            integrity_violations = self.integrity_checker.check_integrity_violations(
                original_protein, optimized_protein
            )
            analysis['protein_integrity'] = integrity_violations
        
        return analysis




