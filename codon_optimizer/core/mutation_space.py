"""
Mutation space for managing possible mutations in DNA sequences.
Inspired by DnaChisel's MutationSpace.
"""

from typing import Dict, List, Tuple, Iterator, Optional
import random
from Bio.Data import CodonTable
from codon_optimizer.core.codon_generator import CodonCombinationGenerator


class MutationSpace:
    """
    Manages the space of possible mutations for a protein sequence.
    
    IMPORTANT: All mutations are SYNONYMOUS - they preserve the amino acid sequence.
    Only codon substitutions are allowed (e.g., AAA -> AAG for Lysine).
    
    Pre-computes all possible codon choices for each amino acid position to enable
    efficient mutation and optimization while maintaining protein integrity.
    
    Example:
        For protein "MK" (Methionine-Lysine):
        - Position 0 (M): Only codon ATG is possible
        - Position 1 (K): Codons AAA, AAG are possible
        - Total variants: 1 * 2 = 2 sequences
    """
    
    def __init__(self, protein_sequence: str, host: str = "CHO"):
        """
        Initialize mutation space from amino acid sequence.
        
        Args:
            protein_sequence: Protein sequence (amino acids, single-letter code)
            host: Host organism (for codon usage preferences)
            
        Raises:
            ValueError: If any amino acid has no codons
        """
        self.protein_sequence = protein_sequence
        self.host = host
        self.codon_generator = CodonCombinationGenerator(host=host)
        self.genetic_code = CodonTable.unambiguous_dna_by_id[1]
        
        # Pre-compute codon options for each position
        self.codon_options: List[List[str]] = []
        for aa in protein_sequence:
            codons = self.codon_generator.get_all_codons_for_aa(aa)
            if not codons:
                raise ValueError(f"No codons found for amino acid: {aa}")
            self.codon_options.append(codons)
        
        # Calculate total space size
        self._space_size = None
    
    @property
    def space_size(self) -> int:
        """Total number of possible sequences."""
        if self._space_size is None:
            total = 1
            for options in self.codon_options:
                total *= len(options)
            self._space_size = total
        return self._space_size
    
    def get_codon_options(self, position: int) -> List[str]:
        """
        Get possible codons for a position.
        
        Args:
            position: Position in protein sequence (0-indexed)
            
        Returns:
            List of possible codons
        """
        if position < 0 or position >= len(self.protein_sequence):
            raise IndexError(f"Position {position} out of range")
        return self.codon_options[position].copy()
    
    def all_variants(self, current_sequence: Optional[str] = None) -> Iterator[str]:
        """
        Generate all possible DNA sequence variants.
        
        IMPORTANT: All variants encode the SAME amino acid sequence (self.protein_sequence).
        Only synonymous codon substitutions are used.
        
        Args:
            current_sequence: Current sequence (ignored, generates all variants)
            
        Yields:
            All possible DNA sequences encoding self.protein_sequence
            
        Example:
            For protein "MK":
            - Yields: "ATGGAA" (ATG=Met, GAA=Lys - wait, GAA is Glu, not Lys!)
            - Actually: "ATGAAA" and "ATGAAG" (both encode Met-Lys)
        """
        import itertools
        for combination in itertools.product(*self.codon_options):
            yield ''.join(combination)
    
    def apply_random_mutations(self, sequence: str, n_mutations: int = 1) -> str:
        """
        Apply random mutations to a sequence.
        
        IMPORTANT: All mutations preserve the amino acid sequence.
        Only synonymous codon substitutions are allowed.
        
        Args:
            sequence: Current DNA sequence (must encode self.protein_sequence)
            n_mutations: Number of codon positions to mutate
            
        Returns:
            Mutated DNA sequence (still encoding the same protein)
        """
        if len(sequence) != len(self.protein_sequence) * 3:
            raise ValueError(
                f"Sequence length ({len(sequence)}) doesn't match protein length "
                f"({len(self.protein_sequence)} * 3 = {len(self.protein_sequence) * 3})"
            )
        
        # Verify sequence encodes correct protein
        from Bio.Seq import Seq
        try:
            translated = str(Seq(sequence).translate())
            if translated != self.protein_sequence:
                raise ValueError(
                    f"Sequence does not encode the expected protein. "
                    f"Expected: {self.protein_sequence[:20]}..., "
                    f"Got: {translated[:20]}..."
                )
        except Exception as e:
            raise ValueError(f"Invalid DNA sequence: {e}")
        
        sequence_list = list(sequence)
        positions = random.sample(range(len(self.protein_sequence)), 
                                 min(n_mutations, len(self.protein_sequence)))
        
        for pos in positions:
            # Get current codon
            codon_start = pos * 3
            current_codon = sequence[codon_start:codon_start + 3]
            
            # Get alternative codons for this amino acid (synonymous mutations only)
            options = self.codon_options[pos]  # All codons encoding self.protein_sequence[pos]
            alternatives = [c for c in options if c != current_codon]
            
            if alternatives:
                new_codon = random.choice(alternatives)
                sequence_list[codon_start:codon_start + 3] = list(new_codon)
        
        mutated_sequence = ''.join(sequence_list)
        
        # Verify mutation preserved protein
        try:
            mutated_protein = str(Seq(mutated_sequence).translate())
            if mutated_protein != self.protein_sequence:
                raise ValueError(
                    f"Mutation changed protein sequence! "
                    f"Expected: {self.protein_sequence[:20]}..., "
                    f"Got: {mutated_protein[:20]}..."
                )
        except Exception as e:
            raise ValueError(f"Mutation verification failed: {e}")
        
        return mutated_sequence
    
    def get_mutations_at_position(self, position: int, current_codon: str) -> List[str]:
        """
        Get all possible mutations at a specific position.
        
        Args:
            position: Position in protein sequence
            current_codon: Current codon at this position
            
        Returns:
            List of alternative codons
        """
        options = self.get_codon_options(position)
        return [c for c in options if c != current_codon]
    
    def create_variant(self, sequence: str, mutations: Dict[int, str]) -> str:
        """
        Create a variant with specific codon mutations.
        
        IMPORTANT: All mutations must preserve the amino acid sequence.
        Only codons from self.codon_options[position] are allowed.
        
        Args:
            sequence: Current DNA sequence
            mutations: Dictionary mapping amino acid position -> new codon
                      (codon must encode the same amino acid as self.protein_sequence[position])
            
        Returns:
            Variant DNA sequence (encoding the same protein)
            
        Raises:
            ValueError: If any mutation would change the amino acid sequence
        """
        sequence_list = list(sequence)
        for pos, new_codon in mutations.items():
            if pos < 0 or pos >= len(self.protein_sequence):
                raise ValueError(f"Position {pos} out of range [0, {len(self.protein_sequence)})")
            if len(new_codon) != 3:
                raise ValueError(f"Invalid codon length: {new_codon} (must be 3 nucleotides)")
            
            # Verify new codon encodes the correct amino acid
            expected_aa = self.protein_sequence[pos]
            if new_codon not in self.codon_options[pos]:
                # Check what amino acid this codon actually encodes
                from Bio.Seq import Seq
                actual_aa = str(Seq(new_codon).translate())
                raise ValueError(
                    f"Invalid codon at position {pos}: '{new_codon}' encodes {actual_aa}, "
                    f"but position {pos} requires {expected_aa}. "
                    f"Valid codons: {self.codon_options[pos]}"
                )
            
            codon_start = pos * 3
            sequence_list[codon_start:codon_start + 3] = list(new_codon)
        
        variant_sequence = ''.join(sequence_list)
        
        # Final verification: translated protein must match
        from Bio.Seq import Seq
        try:
            variant_protein = str(Seq(variant_sequence).translate())
            if variant_protein != self.protein_sequence:
                raise ValueError(
                    f"Variant does not encode the expected protein. "
                    f"Expected: {self.protein_sequence[:30]}..., "
                    f"Got: {variant_protein[:30]}..."
                )
        except Exception as e:
            raise ValueError(f"Variant verification failed: {e}")
        
        return variant_sequence

