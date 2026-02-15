"""
Expression cassette constraints (Kozak, start/stop codons, signal peptides).
"""

from typing import Tuple, Optional, Dict
import re


class ExpressionCassetteChecker:
    """Checker for expression cassette elements."""
    
    def __init__(self,
                 kozak_sequence: str = "GCCACC",
                 start_codon: str = "ATG",
                 stop_codons: Tuple[str, ...] = ("TAA", "TAG", "TGA")):
        """
        Initialize expression cassette checker.
        
        Args:
            kozak_sequence: Kozak consensus sequence
            start_codon: Start codon (usually ATG)
            stop_codons: Valid stop codons
        """
        self.kozak_sequence = kozak_sequence.upper()
        self.start_codon = start_codon.upper()
        self.stop_codons = tuple(c.upper() for c in stop_codons)
    
    def find_start_codon(self, sequence: str) -> Optional[int]:
        """
        Find start codon position.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Position of start codon or None
        """
        sequence = sequence.upper()
        pos = sequence.find(self.start_codon)
        return pos if pos != -1 else None
    
    def find_stop_codon(self, sequence: str, start_pos: int = 0) -> Optional[int]:
        """
        Find stop codon position.
        
        Args:
            sequence: DNA sequence
            start_pos: Position to start searching from
            
        Returns:
            Position of stop codon or None
        """
        sequence = sequence.upper()
        
        for i in range(start_pos, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in self.stop_codons:
                return i
        
        return None
    
    def check_kozak_sequence(self, sequence: str, start_pos: int) -> Dict[str, any]:
        """
        Check Kozak sequence context.
        
        Args:
            sequence: DNA sequence
            start_pos: Position of start codon
            
        Returns:
            Dictionary with Kozak analysis
        """
        if start_pos < len(self.kozak_sequence):
            return {
                'has_kozak': False,
                'kozak_position': None,
                'kozak_sequence': None,
                'kozak_score': 0.0
            }
        
        # Check region before start codon
        kozak_start = start_pos - len(self.kozak_sequence)
        kozak_region = sequence[kozak_start:start_pos].upper()
        
        # Calculate match score
        matches = sum(1 for a, b in zip(kozak_region, self.kozak_sequence) if a == b)
        score = matches / len(self.kozak_sequence)
        
        return {
            'has_kozak': kozak_region == self.kozak_sequence,
            'kozak_position': kozak_start,
            'kozak_sequence': kozak_region,
            'kozak_score': score
        }
    
    def validate_orf(self, sequence: str) -> Dict[str, any]:
        """
        Validate ORF structure (start codon, stop codon, length).
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with validation results
        """
        start_pos = self.find_start_codon(sequence)
        
        if start_pos is None:
            return {
                'is_valid': False,
                'error': 'No start codon found',
                'start_pos': None,
                'stop_pos': None,
                'orf_length': 0
            }
        
        stop_pos = self.find_stop_codon(sequence, start_pos + 3)
        
        if stop_pos is None:
            return {
                'is_valid': False,
                'error': 'No stop codon found',
                'start_pos': start_pos,
                'stop_pos': None,
                'orf_length': 0
            }
        
        orf_length = stop_pos + 3 - start_pos
        
        # Check if length is multiple of 3
        if orf_length % 3 != 0:
            return {
                'is_valid': False,
                'error': f'ORF length ({orf_length}) is not multiple of 3',
                'start_pos': start_pos,
                'stop_pos': stop_pos,
                'orf_length': orf_length
            }
        
        return {
            'is_valid': True,
            'error': None,
            'start_pos': start_pos,
            'stop_pos': stop_pos,
            'orf_length': orf_length
        }
    
    def check_signal_peptide_junction(self,
                                     sequence: str,
                                     signal_peptide_end: int) -> Dict[str, any]:
        """
        Check junction between signal peptide and mature protein.
        
        Args:
            sequence: DNA sequence
            signal_peptide_end: End position of signal peptide
            
        Returns:
            Dictionary with junction analysis
        """
        if signal_peptide_end >= len(sequence):
            return {
                'is_valid': False,
                'error': 'Signal peptide end beyond sequence length'
            }
        
        # Check that junction is at codon boundary
        if signal_peptide_end % 3 != 0:
            return {
                'is_valid': False,
                'error': 'Signal peptide end not at codon boundary',
                'junction_position': signal_peptide_end
            }
        
        # Check for cleavage site patterns (simplified)
        # Real signal peptides have specific cleavage motifs
        junction_codon = sequence[signal_peptide_end:signal_peptide_end + 3]
        
        return {
            'is_valid': True,
            'junction_position': signal_peptide_end,
            'junction_codon': junction_codon
        }
    
    def get_cassette_statistics(self, sequence: str) -> Dict[str, any]:
        """
        Get comprehensive statistics about expression cassette.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with cassette statistics
        """
        orf_validation = self.validate_orf(sequence)
        start_pos = orf_validation.get('start_pos')
        
        kozak_info = {}
        if start_pos is not None:
            kozak_info = self.check_kozak_sequence(sequence, start_pos)
        
        return {
            'orf_valid': orf_validation['is_valid'],
            'start_codon_pos': start_pos,
            'stop_codon_pos': orf_validation.get('stop_pos'),
            'orf_length': orf_validation.get('orf_length', 0),
            'kozak': kozak_info
        }
    
    def calculate_cassette_score(self, sequence: str) -> float:
        """
        Calculate expression cassette score (0-1, higher is better).
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Score between 0 and 1
        """
        orf_validation = self.validate_orf(sequence)
        
        if not orf_validation['is_valid']:
            return 0.0
        
        score = 0.5  # Base score for valid ORF
        
        # Add score for Kozak sequence
        start_pos = orf_validation['start_pos']
        if start_pos is not None:
            kozak_info = self.check_kozak_sequence(sequence, start_pos)
            score += 0.5 * kozak_info.get('kozak_score', 0.0)
        
        return min(1.0, score)




