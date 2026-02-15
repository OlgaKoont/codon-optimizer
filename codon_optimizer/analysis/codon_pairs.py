"""
Codon pair preference analysis for CHO cells.
"""

from typing import Dict, List, Tuple
import numpy as np
from collections import defaultdict


# Codon pair preference matrix for CHO cells
# Based on highly expressed genes - values represent relative frequency
# Higher values indicate preferred pairs
CHO_CODON_PAIR_PREFERENCES = {
    # This is a simplified matrix - in practice, this would be built from
    # analysis of highly expressed CHO genes
    # Format: (codon1, codon2): preference_score
}

# Initialize with default preferences (will be populated from data)
# For now, we'll use a heuristic based on codon usage


class CodonPairAnalyzer:
    """Analyzer for codon pair preferences."""
    
    def __init__(self, host: str = "CHO"):
        """
        Initialize codon pair analyzer.
        
        Args:
            host: Host organism (currently only "CHO" supported)
        """
        self.host = host
        self.pair_matrix = self._build_pair_matrix()
    
    def _build_pair_matrix(self) -> Dict[Tuple[str, str], float]:
        """
        Build codon pair preference matrix.
        In practice, this would be loaded from experimental data.
        For now, we use a heuristic based on individual codon frequencies.
        
        Returns:
            Dictionary mapping (codon1, codon2) tuples to preference scores
        """
        # Import codon usage data
        from codon_optimizer.analysis.codon_usage import CHO_CODON_USAGE
        
        pair_matrix = {}
        
        # Generate all possible codon pairs
        codons = list(CHO_CODON_USAGE.keys())
        
        for codon1 in codons:
            for codon2 in codons:
                # Heuristic: preference based on product of individual codon frequencies
                # In practice, this should be based on observed frequencies in highly expressed genes
                freq1 = CHO_CODON_USAGE.get(codon1, 0.0)
                freq2 = CHO_CODON_USAGE.get(codon2, 0.0)
                
                # Base preference on product, but can be adjusted
                preference = freq1 * freq2
                
                # Add some penalty for certain problematic pairs (e.g., rare tRNA conflicts)
                # This is a placeholder - real data would identify these
                pair_matrix[(codon1, codon2)] = preference
        
        # Normalize to 0-1 range
        max_pref = max(pair_matrix.values()) if pair_matrix else 1.0
        if max_pref > 0:
            pair_matrix = {k: v / max_pref for k, v in pair_matrix.items()}
        
        return pair_matrix
    
    def get_pair_preference(self, codon1: str, codon2: str) -> float:
        """
        Get preference score for a codon pair.
        
        Args:
            codon1: First codon
            codon2: Second codon
            
        Returns:
            Preference score (0-1, higher is better)
        """
        codon1 = codon1.upper()
        codon2 = codon2.upper()
        
        return self.pair_matrix.get((codon1, codon2), 0.5)
    
    def calculate_pair_score(self, dna_sequence: str) -> float:
        """
        Calculate overall codon pair preference score for sequence.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Average pair preference score (0-1)
        """
        dna_sequence = dna_sequence.upper()
        
        if len(dna_sequence) < 6:
            return 0.5  # Not enough codons for pairs
        
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence) - 2, 3)]
        
        if len(codons) < 2:
            return 0.5
        
        pair_scores = []
        for i in range(len(codons) - 1):
            codon1 = codons[i]
            codon2 = codons[i + 1]
            
            if len(codon1) == 3 and len(codon2) == 3:
                preference = self.get_pair_preference(codon1, codon2)
                pair_scores.append(preference)
        
        if not pair_scores:
            return 0.5
        
        return float(np.mean(pair_scores))
    
    def find_poor_pairs(self, dna_sequence: str, threshold: float = 0.3) -> List[Tuple[int, str, str, float]]:
        """
        Find codon pairs with low preference scores.
        
        Args:
            dna_sequence: DNA sequence
            threshold: Minimum acceptable preference score
            
        Returns:
            List of (position, codon1, codon2, score) tuples
        """
        dna_sequence = dna_sequence.upper()
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence) - 2, 3)]
        
        poor_pairs = []
        
        for i in range(len(codons) - 1):
            codon1 = codons[i]
            codon2 = codons[i + 1]
            
            if len(codon1) == 3 and len(codon2) == 3:
                preference = self.get_pair_preference(codon1, codon2)
                
                if preference < threshold:
                    position = i * 3
                    poor_pairs.append((position, codon1, codon2, preference))
        
        return poor_pairs
    
    def suggest_better_pair(self, 
                           codon1: str,
                           codon2: str,
                           amino_acid1: str,
                           amino_acid2: str,
                           codon_usage_analyzer) -> Tuple[str, str]:
        """
        Suggest better codon pair for given amino acids.
        
        Args:
            codon1: Current first codon
            codon2: Current second codon
            amino_acid1: First amino acid
            amino_acid2: Second amino acid
            codon_usage_analyzer: CodonUsageAnalyzer instance
            
        Returns:
            Tuple of (suggested_codon1, suggested_codon2)
        """
        # Get all possible codons for each amino acid
        options1 = codon_usage_analyzer.get_codon_options(amino_acid1)
        options2 = codon_usage_analyzer.get_codon_options(amino_acid2)
        
        # Find best pair
        best_pair = (codon1, codon2)
        best_score = self.get_pair_preference(codon1, codon2)
        
        for opt1 in options1:
            for opt2 in options2:
                score = self.get_pair_preference(opt1, opt2)
                if score > best_score:
                    best_score = score
                    best_pair = (opt1, opt2)
        
        return best_pair
    
    def analyze_context(self, 
                       dna_sequence: str,
                       position: int,
                       context_size: int = 2) -> Dict[str, float]:
        """
        Analyze codon context around a specific position.
        
        Args:
            dna_sequence: DNA sequence
            position: Position to analyze (codon index)
            context_size: Number of codons to consider on each side
            
        Returns:
            Dictionary with context analysis
        """
        dna_sequence = dna_sequence.upper()
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence) - 2, 3)]
        
        if position < 0 or position >= len(codons):
            return {}
        
        context_scores = []
        
        # Analyze pairs before
        for i in range(max(0, position - context_size), position):
            if i + 1 < len(codons):
                codon1 = codons[i]
                codon2 = codons[i + 1]
                if len(codon1) == 3 and len(codon2) == 3:
                    score = self.get_pair_preference(codon1, codon2)
                    context_scores.append(score)
        
        # Analyze pairs after
        for i in range(position, min(len(codons) - 1, position + context_size)):
            codon1 = codons[i]
            codon2 = codons[i + 1]
            if len(codon1) == 3 and len(codon2) == 3:
                score = self.get_pair_preference(codon1, codon2)
                context_scores.append(score)
        
        if not context_scores:
            return {'mean_score': 0.5, 'min_score': 0.5, 'max_score': 0.5}
        
        return {
            'mean_score': float(np.mean(context_scores)),
            'min_score': float(np.min(context_scores)),
            'max_score': float(np.max(context_scores)),
            'std_score': float(np.std(context_scores))
        }
    
    def get_pair_statistics(self, dna_sequence: str) -> Dict[str, float]:
        """
        Get comprehensive statistics about codon pairs in sequence.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Dictionary with pair statistics
        """
        dna_sequence = dna_sequence.upper()
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence) - 2, 3)]
        
        if len(codons) < 2:
            return {
                'mean_preference': 0.5,
                'min_preference': 0.5,
                'max_preference': 0.5,
                'poor_pair_count': 0,
                'poor_pair_fraction': 0.0
            }
        
        pair_scores = []
        for i in range(len(codons) - 1):
            codon1 = codons[i]
            codon2 = codons[i + 1]
            
            if len(codon1) == 3 and len(codon2) == 3:
                preference = self.get_pair_preference(codon1, codon2)
                pair_scores.append(preference)
        
        if not pair_scores:
            return {
                'mean_preference': 0.5,
                'min_preference': 0.5,
                'max_preference': 0.5,
                'poor_pair_count': 0,
                'poor_pair_fraction': 0.0
            }
        
        poor_pairs = [s for s in pair_scores if s < 0.3]
        
        return {
            'mean_preference': float(np.mean(pair_scores)),
            'min_preference': float(np.min(pair_scores)),
            'max_preference': float(np.max(pair_scores)),
            'std_preference': float(np.std(pair_scores)),
            'poor_pair_count': len(poor_pairs),
            'poor_pair_fraction': len(poor_pairs) / len(pair_scores) if pair_scores else 0.0
        }
    
    def load_pair_matrix_from_data(self, pair_frequencies: Dict[Tuple[str, str], int]):
        """
        Load codon pair preference matrix from experimental data.
        
        Args:
            pair_frequencies: Dictionary mapping (codon1, codon2) to observed frequency
        """
        if not pair_frequencies:
            return
        
        # Normalize frequencies
        total = sum(pair_frequencies.values())
        if total > 0:
            self.pair_matrix = {
                pair: count / total for pair, count in pair_frequencies.items()
            }
            
            # Normalize to 0-1 range
            max_freq = max(self.pair_matrix.values()) if self.pair_matrix else 1.0
            if max_freq > 0:
                self.pair_matrix = {k: v / max_freq for k, v in self.pair_matrix.items()}




