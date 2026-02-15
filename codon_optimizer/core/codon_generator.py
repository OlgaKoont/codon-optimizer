"""
Generator for all possible codon combinations for a protein sequence.
"""

from typing import List, Iterator, Tuple
import itertools
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from Bio.Data import CodonTable


class CodonCombinationGenerator:
    """Generator for all possible codon combinations."""
    
    def __init__(self, host: str = "CHO"):
        """
        Initialize codon combination generator.
        
        Args:
            host: Host organism
        """
        self.codon_analyzer = CodonUsageAnalyzer(host=host)
        self.genetic_code = CodonTable.unambiguous_dna_by_id[1]
    
    def get_all_codons_for_aa(self, amino_acid: str) -> List[str]:
        """
        Get all possible codons for an amino acid.
        
        Args:
            amino_acid: Single letter amino acid code
            
        Returns:
            List of all codons encoding this amino acid
        """
        codons = []
        for codon, aa in self.genetic_code.forward_table.items():
            if aa == amino_acid:
                codons.append(codon)
        return codons
    
    def generate_all_combinations(self, protein_sequence: str) -> Iterator[str]:
        """
        Generate all possible DNA sequences for a protein sequence.
        
        WARNING: This can generate a huge number of combinations!
        For a protein of length N with average M codons per amino acid,
        total combinations = M^N
        
        Args:
            protein_sequence: Protein sequence (amino acids)
            
        Yields:
            DNA sequences encoding the protein
        """
        # Get codon options for each amino acid
        codon_options_list = []
        for aa in protein_sequence:
            codons = self.get_all_codons_for_aa(aa)
            if not codons:
                raise ValueError(f"No codons found for amino acid: {aa}")
            codon_options_list.append(codons)
        
        # Generate all combinations using itertools.product
        for combination in itertools.product(*codon_options_list):
            yield ''.join(combination)
    
    def estimate_combination_count(self, protein_sequence: str) -> int:
        """
        Estimate total number of possible combinations.
        
        Args:
            protein_sequence: Protein sequence
            
        Returns:
            Estimated number of combinations
        """
        total = 1
        for aa in protein_sequence:
            codons = self.get_all_codons_for_aa(aa)
            total *= len(codons)
        return total
    
    def generate_smart_combinations(self,
                                   protein_sequence: str,
                                   max_combinations: int = 10000,
                                   strategy: str = 'diverse') -> Iterator[str]:
        """
        Generate a smart subset of codon combinations.
        
        For long sequences, generating all combinations is infeasible.
        This method generates a diverse subset.
        
        Args:
            protein_sequence: Protein sequence
            max_combinations: Maximum number of combinations to generate
            strategy: Strategy for selection ('diverse', 'optimal', 'random')
            
        Yields:
            DNA sequences encoding the protein
        """
        # Get codon options for each amino acid
        codon_options_list = []
        for aa in protein_sequence:
            codons = self.get_all_codons_for_aa(aa)
            if not codons:
                raise ValueError(f"No codons found for amino acid: {aa}")
            codon_options_list.append(codons)
        
        total_combinations = self.estimate_combination_count(protein_sequence)
        
        if total_combinations <= max_combinations:
            # Generate all if feasible
            yield from self.generate_all_combinations(protein_sequence)
        else:
            # Generate diverse subset
            if strategy == 'optimal':
                # Use optimal codons with some variation
                yield from self._generate_optimal_variants(
                    protein_sequence, codon_options_list, max_combinations
                )
            elif strategy == 'diverse':
                # Generate diverse combinations
                yield from self._generate_diverse_combinations(
                    codon_options_list, max_combinations
                )
            else:  # random
                # Random sampling
                yield from self._generate_random_combinations(
                    codon_options_list, max_combinations
                )
    
    def _generate_optimal_variants(self,
                                  protein_sequence: str,
                                  codon_options_list: List[List[str]],
                                  max_combinations: int) -> Iterator[str]:
        """Generate variants starting from optimal codons."""
        import random
        
        # Start with optimal sequence
        optimal_sequence = []
        for i, aa in enumerate(protein_sequence):
            optimal_codon = self.codon_analyzer.get_optimal_codon(aa)
            optimal_sequence.append(optimal_codon)
        
        yield ''.join(optimal_sequence)
        
        # Generate variants by mutating positions
        generated = set([''.join(optimal_sequence)])
        
        for _ in range(max_combinations - 1):
            variant = optimal_sequence.copy()
            
            # Mutate random positions
            num_mutations = min(3, len(protein_sequence))
            positions = random.sample(range(len(protein_sequence)), num_mutations)
            
            for pos in positions:
                aa = protein_sequence[pos]
                options = codon_options_list[pos]
                # Choose different codon
                new_codon = random.choice([c for c in options if c != variant[pos]])
                variant[pos] = new_codon
            
            variant_str = ''.join(variant)
            if variant_str not in generated:
                generated.add(variant_str)
                yield variant_str
    
    def _generate_diverse_combinations(self,
                                      codon_options_list: List[List[str]],
                                      max_combinations: int) -> Iterator[str]:
        """Generate diverse combinations using systematic sampling."""
        import random
        
        generated = set()
        
        # Try different strategies to get diversity
        strategies = ['optimal', 'balanced', 'random']
        
        for strategy in strategies:
            for _ in range(max_combinations // len(strategies)):
                combination = []
                for options in codon_options_list:
                    if strategy == 'optimal':
                        # Prefer first (usually optimal) codon
                        codon = options[0] if options else 'NNN'
                    elif strategy == 'balanced':
                        # Middle codon
                        codon = options[len(options) // 2] if options else 'NNN'
                    else:  # random
                        codon = random.choice(options) if options else 'NNN'
                    combination.append(codon)
                
                combo_str = ''.join(combination)
                if combo_str not in generated:
                    generated.add(combo_str)
                    yield combo_str
                    
                    if len(generated) >= max_combinations:
                        return
    
    def _generate_random_combinations(self,
                                     codon_options_list: List[List[str]],
                                     max_combinations: int) -> Iterator[str]:
        """Generate random combinations."""
        import random
        
        generated = set()
        
        for _ in range(max_combinations):
            combination = []
            for options in codon_options_list:
                codon = random.choice(options) if options else 'NNN'
                combination.append(codon)
            
            combo_str = ''.join(combination)
            if combo_str not in generated:
                generated.add(combo_str)
                yield combo_str


