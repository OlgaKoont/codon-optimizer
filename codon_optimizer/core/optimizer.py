"""
Optimizer for codon optimization from protein sequences.
Generates all possible codon combinations and selects best based on criteria.
"""

import random
from typing import List, Tuple, Dict, Optional
from Bio.Seq import Seq
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.core.scorer import MultiCriteriaScorer
from codon_optimizer.core.sequence import OptimizedSequence
from codon_optimizer.core.codon_generator import CodonCombinationGenerator
from codon_optimizer.core.local_optimizer import LocalCodonOptimizer
from codon_optimizer.core.optimizer_v2 import CodonOptimizerV2
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer


class CodonOptimizer:
    """Optimizer for codon optimization from protein sequences."""
    
    def __init__(self, config: OptimizationConfig):
        """
        Initialize optimizer.
        
        Args:
            config: Optimization configuration
        """
        self.config = config
        self.scorer = MultiCriteriaScorer(config)
        self.codon_analyzer = CodonUsageAnalyzer(
            host=config.host,
            max_codon_usage=config.max_codon_usage
        )
        self.codon_generator = CodonCombinationGenerator(host=config.host)
        self.local_optimizer = LocalCodonOptimizer(config)
        self.optimizer_v2 = CodonOptimizerV2(config)  # New architecture
    
    def _dna_to_protein(self, dna_sequence: str) -> str:
        """Convert DNA to protein sequence."""
        return str(Seq(dna_sequence).translate())
    
    def _generate_initial_population(self,
                                     original_dna: str,
                                     population_size: int) -> List[str]:
        """
        Generate initial population of codon-optimized sequences.
        
        Args:
            original_dna: Original DNA sequence
            population_size: Size of population
            
        Returns:
            List of DNA sequences
        """
        original_protein = self._dna_to_protein(original_dna)
        population = []
        
        for _ in range(population_size):
            optimized = self._optimize_codons_for_protein(
                original_protein,
                strategy=random.choice(['random', 'optimal', 'balanced'])
            )
            population.append(optimized)
        
        return population
    
    def _optimize_codons_for_protein(self,
                                    protein_sequence: str,
                                    strategy: str = 'balanced') -> str:
        """
        Optimize codons for a protein sequence using a strategy.
        
        Args:
            protein_sequence: Protein sequence
            strategy: Optimization strategy ('random', 'optimal', 'balanced')
            
        Returns:
            Optimized DNA sequence
        """
        codons = []
        
        for aa in protein_sequence:
            if strategy == 'random':
                # Random codon choice
                options = self.codon_analyzer.get_codon_options(aa)
                codon = random.choice(options) if options else 'NNN'
            elif strategy == 'optimal':
                # Always use optimal codon
                codon = self.codon_analyzer.get_optimal_codon(aa)
            else:  # balanced
                # Prefer optimal but with some diversity
                options = self.codon_analyzer.get_codon_options(aa, avoid_overuse=True)
                # Weight towards optimal but allow alternatives
                if len(options) > 1:
                    weights = [0.5 if i == 0 else 0.5 / (len(options) - 1) 
                              for i in range(len(options))]
                    codon = random.choices(options, weights=weights)[0]
                else:
                    codon = options[0] if options else 'NNN'
            
            codons.append(codon)
        
        return ''.join(codons)
    
    def _mutate(self,
               dna_sequence: str,
               mutation_rate: float) -> str:
        """
        Mutate sequence by changing codons while preserving protein.
        
        Args:
            dna_sequence: DNA sequence
            mutation_rate: Probability of mutating each codon
            
        Returns:
            Mutated DNA sequence
        """
        protein = self._dna_to_protein(dna_sequence)
        codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
        
        mutated_codons = []
        for i, (codon, aa) in enumerate(zip(codons, protein)):
            if random.random() < mutation_rate:
                # Mutate this codon
                options = self.codon_analyzer.get_codon_options(aa, avoid_overuse=True)
                # Prefer alternatives to current codon
                alternatives = [c for c in options if c != codon]
                if alternatives:
                    new_codon = random.choice(alternatives)
                else:
                    new_codon = codon
                mutated_codons.append(new_codon)
            else:
                mutated_codons.append(codon)
        
        return ''.join(mutated_codons)
    
    def _crossover(self,
                  parent1: str,
                  parent2: str) -> Tuple[str, str]:
        """
        Perform crossover between two parent sequences.
        
        Args:
            parent1: First parent DNA sequence
            parent2: Second parent DNA sequence
            
        Returns:
            Tuple of (offspring1, offspring2)
        """
        # Ensure same protein sequence
        protein1 = self._dna_to_protein(parent1)
        protein2 = self._dna_to_protein(parent2)
        
        if protein1 != protein2:
            return parent1, parent2
        
        # Crossover at codon boundaries
        codons1 = [parent1[i:i+3] for i in range(0, len(parent1), 3)]
        codons2 = [parent2[i:i+3] for i in range(0, len(parent2), 3)]
        
        crossover_point = random.randint(1, len(codons1) - 1)
        
        offspring1_codons = codons1[:crossover_point] + codons2[crossover_point:]
        offspring2_codons = codons2[:crossover_point] + codons1[crossover_point:]
        
        return ''.join(offspring1_codons), ''.join(offspring2_codons)
    
    def _select_parents(self,
                       population: List[str],
                       scores: List[float],
                       num_parents: int) -> List[str]:
        """
        Select parents using tournament selection.
        
        Args:
            population: List of DNA sequences
            scores: List of fitness scores
            num_parents: Number of parents to select
            
        Returns:
            List of selected parent sequences
        """
        parents = []
        tournament_size = 3
        
        for _ in range(num_parents):
            # Tournament selection
            tournament_indices = random.sample(range(len(population)), 
                                             min(tournament_size, len(population)))
            tournament_scores = [(i, scores[i]) for i in tournament_indices]
            tournament_scores.sort(key=lambda x: x[1], reverse=True)
            winner_idx = tournament_scores[0][0]
            parents.append(population[winner_idx])
        
        return parents
    
    def optimize(self,
                protein_sequence: str,
                max_combinations: Optional[int] = None,
                use_local_optimization: Optional[bool] = None,
                progress_callback=None) -> OptimizedSequence:
        """
        Optimize codon usage from protein sequence.
        
        Uses two strategies:
        1. Full enumeration for short sequences (all combinations)
        2. Local iterative optimization for long sequences (DnaChisel-inspired)
        
        Args:
            protein_sequence: Protein sequence (amino acids)
            max_combinations: Maximum number of combinations to evaluate
                             (None for config default or all if feasible)
            use_local_optimization: Force use of local optimization (None = auto-detect)
            
        Returns:
            Optimized sequence
        """
        # Estimate total combinations
        total_combinations = self.codon_generator.estimate_combination_count(protein_sequence)
        
        # Try new architecture first (DnaChisel-inspired)
        try:
            return self.optimizer_v2.optimize(protein_sequence, max_combinations, progress_callback=progress_callback)
        except Exception:
            # Fallback to old methods if new architecture fails
            pass
        
        # Auto-detect strategy if not specified
        if use_local_optimization is None:
            # Use local optimization for sequences longer than 100 aa or if combinations > 10000
            use_local_optimization = (
                len(protein_sequence) > 100 or 
                total_combinations > 10000
            )
        
        if use_local_optimization:
            # Use local iterative optimization (DnaChisel approach)
            optimized_dna, metadata = self.local_optimizer.optimize(protein_sequence)
            
            # Calculate final score
            score_dict = self.scorer.calculate_total_score(optimized_dna, protein_sequence)
            
            return OptimizedSequence(
                dna_sequence=optimized_dna,
                original_sequence=None,
                sequence_id="optimized",
                metadata={
                    'fitness_score': score_dict['total_score'],
                    'optimization_method': 'local_iterative',
                    **metadata
                }
            )
        else:
            # Use full enumeration (original approach)
            # Determine max combinations
            if max_combinations is None:
                max_combinations = min(total_combinations, self.config.population_size * 100)
            
            # Generate codon combinations
            if total_combinations <= max_combinations:
                # Generate all combinations if feasible
                combinations = list(self.codon_generator.generate_all_combinations(protein_sequence))
            else:
                # Generate smart subset
                combinations = list(self.codon_generator.generate_smart_combinations(
                    protein_sequence,
                    max_combinations=max_combinations,
                    strategy='diverse'
                ))
            
            # Evaluate all combinations
            best_sequence = None
            best_score = -float('inf')
            scored_sequences = []
            
            for dna_seq in combinations:
                score_dict = self.scorer.calculate_total_score(
                    dna_seq, protein_sequence
                )
                score = score_dict['total_score']
                
                scored_sequences.append((dna_seq, score, score_dict))
                
                if score > best_score:
                    best_score = score
                    best_sequence = dna_seq
            
            if best_sequence is None:
                raise ValueError("No valid codon combinations found")
            
            return OptimizedSequence(
                dna_sequence=best_sequence,
                original_sequence=None,  # No original DNA, only protein
                sequence_id="optimized",
                metadata={
                    'fitness_score': best_score,
                    'total_combinations_evaluated': len(combinations),
                    'total_possible_combinations': total_combinations,
                    'optimization_method': 'full_enumeration'
                }
            )
    
    def optimize_multiple(self,
                         protein_sequence: str,
                         num_solutions: int = 10,
                         max_combinations: Optional[int] = None) -> List[OptimizedSequence]:
        """
        Generate multiple optimized solutions from protein sequence.
        
        Args:
            protein_sequence: Protein sequence (amino acids)
            num_solutions: Number of top solutions to return
            max_combinations: Maximum number of combinations to evaluate
            
        Returns:
            List of optimized sequences, sorted by score
        """
        # Estimate total combinations
        total_combinations = self.codon_generator.estimate_combination_count(protein_sequence)
        
        # Determine max combinations
        if max_combinations is None:
            max_combinations = min(total_combinations, self.config.population_size * 100)
        
        # Generate codon combinations
        if total_combinations <= max_combinations:
            combinations = list(self.codon_generator.generate_all_combinations(protein_sequence))
        else:
            combinations = list(self.codon_generator.generate_smart_combinations(
                protein_sequence,
                max_combinations=max_combinations,
                strategy='diverse'
            ))
        
        # Evaluate all combinations
        scored_sequences = []
        
        for dna_seq in combinations:
            score_dict = self.scorer.calculate_total_score(
                dna_seq, protein_sequence
            )
            score = score_dict['total_score']
            scored_sequences.append((dna_seq, score, score_dict))
        
        # Sort by score and take top N
        scored_sequences.sort(key=lambda x: x[1], reverse=True)
        
        solutions = []
        for i, (dna_seq, score, score_dict) in enumerate(scored_sequences[:num_solutions]):
            solution = OptimizedSequence(
                dna_sequence=dna_seq,
                original_sequence=None,
                sequence_id=f"solution_{i+1}",
                metadata={
                    'fitness_score': score,
                    'rank': i + 1,
                    'total_combinations_evaluated': len(combinations),
                    'total_possible_combinations': total_combinations
                }
            )
            solutions.append(solution)
        
        return solutions



