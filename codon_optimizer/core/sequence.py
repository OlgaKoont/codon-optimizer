"""
Base sequence class for codon optimization.
"""

from typing import List, Dict, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import seq1


class OptimizedSequence:
    """Represents an optimized DNA sequence with metadata."""
    
    def __init__(self, 
                 dna_sequence: str,
                 original_sequence: Optional[str] = None,
                 original_protein: Optional[str] = None,
                 sequence_id: str = "optimized",
                 metadata: Optional[Dict] = None):
        """
        Initialize optimized sequence.
        
        Args:
            dna_sequence: Optimized DNA sequence
            original_sequence: Original DNA sequence (optional)
            original_protein: Original protein sequence (optional, for validation)
            sequence_id: Sequence identifier
            metadata: Additional metadata dictionary
        """
        self.dna_sequence = dna_sequence.upper()
        self.original_sequence = original_sequence.upper() if original_sequence else None
        self.original_protein = original_protein
        self.sequence_id = sequence_id
        self.metadata = metadata or {}
        
        # Validate sequences have same protein sequence if original provided
        if self.original_sequence:
            self._validate_protein_integrity()
        elif self.original_protein:
            self._validate_protein_integrity_from_protein()
    
    def _validate_protein_integrity(self):
        """Validate that optimized sequence encodes same protein as original DNA."""
        original_protein = str(Seq(self.original_sequence).translate())
        optimized_protein = str(Seq(self.dna_sequence).translate())
        
        if original_protein != optimized_protein:
            raise ValueError(
                f"Protein sequences differ!\n"
                f"Original: {original_protein}\n"
                f"Optimized: {optimized_protein}"
            )
    
    def _validate_protein_integrity_from_protein(self):
        """Validate that optimized sequence encodes same protein as original protein."""
        optimized_protein = str(Seq(self.dna_sequence).translate())
        
        if optimized_protein != self.original_protein:
            raise ValueError(
                f"Protein sequences differ!\n"
                f"Original: {self.original_protein}\n"
                f"Optimized: {optimized_protein}"
            )
    
    def get_protein_sequence(self) -> str:
        """Get protein sequence from DNA sequence."""
        return str(Seq(self.dna_sequence).translate())
    
    def get_codons(self) -> List[str]:
        """Get list of codons from DNA sequence."""
        codons = []
        for i in range(0, len(self.dna_sequence) - 2, 3):
            codon = self.dna_sequence[i:i+3]
            if len(codon) == 3:
                codons.append(codon)
        return codons
    
    def get_amino_acids(self) -> List[str]:
        """Get list of amino acids from DNA sequence."""
        codons = self.get_codons()
        return [str(Seq(codon).translate()) for codon in codons]
    
    def __len__(self) -> int:
        """Return length of DNA sequence."""
        return len(self.dna_sequence)
    
    def __str__(self) -> str:
        """Return DNA sequence as string."""
        return self.dna_sequence
    
    def __repr__(self) -> str:
        """Return representation of sequence."""
        return f"OptimizedSequence(id={self.sequence_id}, length={len(self)})"


class SequencePair:
    """Represents a pair of sequences (e.g., heavy and light chains)."""
    
    def __init__(self,
                 sequence1: OptimizedSequence,
                 sequence2: OptimizedSequence,
                 pair_type: str = "HC_LC"):
        """
        Initialize sequence pair.
        
        Args:
            sequence1: First sequence (e.g., heavy chain)
            sequence2: Second sequence (e.g., light chain)
            pair_type: Type of pair (e.g., "HC_LC")
        """
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.pair_type = pair_type
    
    def get_total_length(self) -> int:
        """Get total length of both sequences."""
        return len(self.sequence1) + len(self.sequence2)
    
    def __repr__(self) -> str:
        """Return representation of pair."""
        return f"SequencePair(type={self.pair_type}, len1={len(self.sequence1)}, len2={len(self.sequence2)})"



