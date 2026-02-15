"""
Codon usage analysis and CAI calculation for CHO cells.
"""

from typing import Dict, List, Tuple
import numpy as np
from Bio.Data import CodonTable


# CHO codon usage table based on highly expressed genes
# Values represent relative frequency of each codon (normalized to 1.0 for most frequent codon per amino acid)
CHO_CODON_USAGE = {
    'TTT': 0.58, 'TTC': 1.00,  # Phe
    'TTA': 0.12, 'TTG': 0.15, 'CTT': 0.18, 'CTC': 0.35, 'CTA': 0.08, 'CTG': 1.00,  # Leu
    'ATT': 0.25, 'ATC': 0.75, 'ATA': 0.08,  # Ile
    'ATG': 1.00,  # Met
    'GTT': 0.28, 'GTC': 0.45, 'GTA': 0.12, 'GTG': 1.00,  # Val
    'TCT': 0.18, 'TCC': 0.35, 'TCA': 0.15, 'TCG': 0.08, 'AGT': 0.12, 'AGC': 0.45,  # Ser
    'CCT': 0.18, 'CCC': 0.35, 'CCA': 0.45, 'CCG': 0.12,  # Pro
    'ACT': 0.18, 'ACC': 0.75, 'ACA': 0.25, 'ACG': 0.12,  # Thr
    'GCT': 0.25, 'GCC': 0.75, 'GCA': 0.18, 'GCG': 0.12,  # Ala
    'TAT': 0.45, 'TAC': 1.00,  # Tyr
    'CAT': 0.35, 'CAC': 1.00,  # His
    'CAA': 0.45, 'CAG': 1.00,  # Gln
    'AAT': 0.35, 'AAC': 1.00,  # Asn
    'AAA': 0.45, 'AAG': 1.00,  # Lys
    'GAT': 0.45, 'GAC': 1.00,  # Asp
    'GAA': 0.65, 'GAG': 1.00,  # Glu
    'TGT': 0.35, 'TGC': 1.00,  # Cys
    'TGG': 1.00,  # Trp
    'CGT': 0.12, 'CGC': 0.35, 'CGA': 0.08, 'CGG': 0.15, 'AGA': 0.18, 'AGG': 0.45,  # Arg
    'GGT': 0.25, 'GGC': 0.75, 'GGA': 0.18, 'GGG': 0.35,  # Gly
    'TAA': 1.00, 'TAG': 0.35, 'TGA': 0.65,  # Stop
}

# Standard genetic code
GENETIC_CODE = CodonTable.unambiguous_dna_by_id[1]


class CodonUsageAnalyzer:
    """Analyzer for codon usage and CAI calculation."""
    
    def __init__(self, host: str = "CHO", max_codon_usage: float = 0.5):
        """
        Initialize codon usage analyzer.
        
        Args:
            host: Host organism (currently only "CHO" supported)
            max_codon_usage: Maximum allowed usage of single codon (to avoid tRNA depletion)
        """
        self.host = host
        self.max_codon_usage = max_codon_usage
        
        if host == "CHO":
            self.codon_table = CHO_CODON_USAGE.copy()
        else:
            raise ValueError(f"Host {host} not yet supported. Only 'CHO' is available.")
        
        # Build amino acid to codons mapping
        self.aa_to_codons = self._build_aa_to_codons()
        
        # Calculate relative adaptiveness (w values) for CAI
        self.w_values = self._calculate_w_values()
    
    def _build_aa_to_codons(self) -> Dict[str, List[str]]:
        """Build mapping from amino acid to list of codons."""
        aa_to_codons = {}
        for codon, aa in GENETIC_CODE.forward_table.items():
            if aa not in aa_to_codons:
                aa_to_codons[aa] = []
            aa_to_codons[aa].append(codon)
        return aa_to_codons
    
    def _calculate_w_values(self) -> Dict[str, float]:
        """
        Calculate relative adaptiveness (w) for each codon.
        w = frequency of codon / frequency of most frequent codon for that amino acid
        """
        w_values = {}
        
        for aa, codons in self.aa_to_codons.items():
            # Find maximum frequency for this amino acid
            max_freq = max(self.codon_table.get(codon, 0.0) for codon in codons)
            
            if max_freq > 0:
                for codon in codons:
                    freq = self.codon_table.get(codon, 0.0)
                    w_values[codon] = freq / max_freq if max_freq > 0 else 0.0
            else:
                # If no frequency data, set all to 1.0
                for codon in codons:
                    w_values[codon] = 1.0
        
        return w_values
    
    def calculate_cai(self, dna_sequence: str) -> float:
        """
        Calculate Codon Adaptation Index (CAI) for a DNA sequence.
        
        CAI = (w1 * w2 * ... * wn)^(1/n)
        where wi is the relative adaptiveness of codon i
        
        Args:
            dna_sequence: DNA sequence (must be multiple of 3)
            
        Returns:
            CAI value between 0 and 1
        """
        dna_sequence = dna_sequence.upper()
        
        if len(dna_sequence) % 3 != 0:
            raise ValueError("Sequence length must be multiple of 3")
        
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
        
        # Filter out stop codons
        codons = [c for c in codons if c not in ['TAA', 'TAG', 'TGA']]
        
        if not codons:
            return 0.0
        
        # Calculate product of w values
        w_product = 1.0
        for codon in codons:
            w = self.w_values.get(codon, 0.0)
            if w == 0.0:
                # Unknown codon, use minimum value
                w = 0.01
            w_product *= w
        
        # Calculate geometric mean
        cai = w_product ** (1.0 / len(codons))
        
        return cai
    
    def get_codon_frequencies(self, dna_sequence: str) -> Dict[str, int]:
        """
        Get codon frequency counts from sequence.
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Dictionary mapping codons to counts
        """
        dna_sequence = dna_sequence.upper()
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence) - 2, 3)]
        
        freq = {}
        for codon in codons:
            if len(codon) == 3:
                freq[codon] = freq.get(codon, 0) + 1
        
        return freq
    
    def check_tRNA_depletion(self, dna_sequence: str) -> Tuple[bool, Dict[str, float]]:
        """
        Check if sequence uses any codon too frequently (risk of tRNA depletion).
        
        Args:
            dna_sequence: DNA sequence
            
        Returns:
            Tuple of (has_depletion_risk, codon_usage_ratios)
        """
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence) - 2, 3)]
        total_codons = len(codons)
        
        if total_codons == 0:
            return False, {}
        
        codon_counts = {}
        for codon in codons:
            if len(codon) == 3:
                codon_counts[codon] = codon_counts.get(codon, 0) + 1
        
        usage_ratios = {codon: count / total_codons for codon, count in codon_counts.items()}
        
        has_risk = any(ratio > self.max_codon_usage for ratio in usage_ratios.values())
        
        return has_risk, usage_ratios
    
    def get_optimal_codon(self, amino_acid: str) -> str:
        """
        Get the most frequently used codon for an amino acid.
        
        Args:
            amino_acid: Single letter amino acid code
            
        Returns:
            Most optimal codon
        """
        codons = self.aa_to_codons.get(amino_acid, [])
        if not codons:
            raise ValueError(f"No codons found for amino acid: {amino_acid}")
        
        # Find codon with highest frequency
        best_codon = codons[0]
        best_freq = self.codon_table.get(best_codon, 0.0)
        
        for codon in codons[1:]:
            freq = self.codon_table.get(codon, 0.0)
            if freq > best_freq:
                best_freq = freq
                best_codon = codon
        
        return best_codon
    
    def get_codon_options(self, amino_acid: str, avoid_overuse: bool = True) -> List[str]:
        """
        Get list of codons for an amino acid, optionally avoiding overused ones.
        
        Args:
            amino_acid: Single letter amino acid code
            avoid_overuse: If True, prefer codons that are less likely to be overused
            
        Returns:
            List of codons sorted by preference
        """
        codons = self.aa_to_codons.get(amino_acid, [])
        if not codons:
            return []
        
        # Sort by frequency (descending)
        codons_with_freq = [(codon, self.codon_table.get(codon, 0.0)) for codon in codons]
        codons_with_freq.sort(key=lambda x: x[1], reverse=True)
        
        if avoid_overuse:
            # Prefer codons with moderate frequency to avoid depletion
            # But still prioritize high-frequency codons
            codons_with_freq.sort(key=lambda x: (x[1] > 0.5, x[1]), reverse=True)
        
        return [codon for codon, _ in codons_with_freq]
    
    def suggest_codon_replacement(self, 
                                  current_codon: str, 
                                  current_usage_ratio: float) -> str:
        """
        Suggest alternative codon if current one is overused.
        
        Args:
            current_codon: Current codon
            current_usage_ratio: Current usage ratio in sequence
            
        Returns:
            Suggested alternative codon
        """
        if current_usage_ratio <= self.max_codon_usage:
            return current_codon
        
        # Get amino acid for this codon
        amino_acid = GENETIC_CODE.forward_table.get(current_codon)
        if not amino_acid:
            return current_codon
        
        # Get alternative codons
        alternatives = self.get_codon_options(amino_acid, avoid_overuse=True)
        
        # Prefer alternatives that are less used but still efficient
        for alt_codon in alternatives:
            if alt_codon != current_codon:
                alt_freq = self.codon_table.get(alt_codon, 0.0)
                # Prefer codons with good frequency but not the most common
                if 0.3 <= alt_freq <= 0.8:
                    return alt_codon
        
        # If no good alternative, return second best
        if len(alternatives) > 1:
            return alternatives[1] if alternatives[0] == current_codon else alternatives[0]
        
        return current_codon



