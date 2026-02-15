"""
Pareto optimization for multi-objective codon optimization.
"""

from typing import List, Dict, Tuple
import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.core.scorer import MultiCriteriaScorer
from codon_optimizer.core.sequence import OptimizedSequence
from codon_optimizer.core.optimizer import CodonOptimizer
from Bio.Seq import Seq


class CodonOptimizationProblem(Problem):
    """Pareto optimization problem for codon optimization."""
    
    def __init__(self,
                 original_protein: str,
                 scorer: MultiCriteriaScorer,
                 codon_analyzer,
                 num_codons: int):
        """
        Initialize optimization problem.
        
        Args:
            original_protein: Original protein sequence (amino acids)
            scorer: MultiCriteriaScorer instance
            codon_analyzer: CodonUsageAnalyzer instance
            num_codons: Number of codons in sequence
        """
        self.original_protein = original_protein
        self.scorer = scorer
        self.codon_analyzer = codon_analyzer
        
        # Build codon options for each amino acid
        self.codon_options = []
        for aa in original_protein:
            options = codon_analyzer.get_codon_options(aa, avoid_overuse=True)
            self.codon_options.append(options)
        
        # Problem definition: minimize two objectives
        # Objective 1: Negative CAI (minimize = maximize CAI)
        # Objective 2: Motif penalty (minimize)
        n_var = num_codons  # Decision variables: codon indices
        n_obj = 2  # Two objectives
        n_constr = 0  # No constraints (protein integrity handled in evaluation)
        
        # Bounds: each variable is index into codon options
        xl = np.zeros(n_var, dtype=int)
        xu = np.array([len(options) - 1 for options in self.codon_options], dtype=int)
        
        super().__init__(n_var=n_var, n_obj=n_obj, n_constr=n_constr, xl=xl, xu=xu, type_var=int)
    
    def _decode_solution(self, x: np.ndarray) -> str:
        """
        Decode solution vector to DNA sequence.
        
        Args:
            x: Solution vector (codon indices)
            
        Returns:
            DNA sequence
        """
        codons = []
        for i, codon_idx in enumerate(x):
            options = self.codon_options[i]
            codon = options[codon_idx] if codon_idx < len(options) else options[0]
            codons.append(codon)
        
        return ''.join(codons)
    
    def _evaluate(self, x, out, *args, **kwargs):
        """
        Evaluate solutions.
        
        Args:
            x: Solution matrix (population x variables)
            out: Output dictionary
        """
        n_pop = x.shape[0]
        objectives = np.zeros((n_pop, 2))
        
        for i in range(n_pop):
            dna_sequence = self._decode_solution(x[i])
            
            # Validate protein integrity
            try:
                optimized_protein = str(Seq(dna_sequence).translate())
                if optimized_protein != self.original_protein:
                    # Invalid solution - high penalty
                    objectives[i, 0] = 1.0  # Bad CAI
                    objectives[i, 1] = 1.0  # Bad motifs
                    continue
            except:
                objectives[i, 0] = 1.0
                objectives[i, 1] = 1.0
                continue
            
            # Calculate objectives
            scores = self.scorer.calculate_total_score(
                dna_sequence, self.original_protein
            )
            
            # Objective 1: Negative CAI (minimize = maximize CAI)
            cai = self.codon_analyzer.calculate_cai(dna_sequence)
            objectives[i, 0] = 1.0 - cai  # Minimize this
            
            # Objective 2: Motif penalty
            motif_penalty = self.scorer.calculate_motif_score(dna_sequence)
            objectives[i, 1] = 1.0 - motif_penalty  # Minimize this
        
        out["F"] = objectives


class ParetoOptimizer:
    """Pareto optimizer for multi-objective codon optimization."""
    
    def __init__(self, config: OptimizationConfig):
        """
        Initialize Pareto optimizer.
        
        Args:
            config: Optimization configuration
        """
        self.config = config
        self.scorer = MultiCriteriaScorer(config)
        from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
        self.codon_analyzer = CodonUsageAnalyzer(
            host=config.host,
            max_codon_usage=config.max_codon_usage
        )
    
    def optimize_pareto(self,
                       protein_sequence: str,
                       num_solutions: int = 10) -> List[OptimizedSequence]:
        """
        Generate Pareto-optimal solutions from protein sequence.
        
        Args:
            protein_sequence: Protein sequence (amino acids)
            num_solutions: Number of solutions to generate
            
        Returns:
            List of Pareto-optimal sequences
        """
        num_codons = len(protein_sequence)
        
        # Create problem
        problem = CodonOptimizationProblem(
            original_protein=protein_sequence,
            scorer=self.scorer,
            codon_analyzer=self.codon_analyzer,
            num_codons=num_codons
        )
        
        # Initialize algorithm
        algorithm = NSGA2(pop_size=num_solutions * 2)
        
        # Optimize
        res = minimize(
            problem,
            algorithm,
            ('n_gen', 50),
            verbose=False,
            seed=42
        )
        
        # Extract Pareto front
        pareto_solutions = []
        
        for i, x in enumerate(res.X):
            dna_sequence = problem._decode_solution(x)
            
            # Validate
            try:
                optimized_protein = str(Seq(dna_sequence).translate())
                if optimized_protein != original_protein:
                    continue
            except:
                continue
            
            # Calculate full scores
            scores = self.scorer.calculate_total_score(dna_sequence, original_protein)
            
            solution = OptimizedSequence(
                dna_sequence=dna_sequence,
                original_sequence=None,  # No original DNA, only protein
                sequence_id=f"pareto_{i+1}",
                metadata={
                    'fitness_score': scores['total_score'],
                    'cai': self.codon_analyzer.calculate_cai(dna_sequence),
                    'motif_score': scores['motif_score'],
                    'pareto_rank': i + 1
                }
            )
            
            pareto_solutions.append(solution)
        
        # Sort by total score
        pareto_solutions.sort(
            key=lambda s: s.metadata.get('fitness_score', 0.0),
            reverse=True
        )
        
        return pareto_solutions[:num_solutions]
    
    def optimize_expression_vs_structure(self,
                                        protein_sequence: str,
                                        num_solutions: int = 10) -> List[OptimizedSequence]:
        """
        Optimize for trade-off between expression (CAI) and structure (mRNA stability).
        
        Args:
            protein_sequence: Protein sequence (amino acids)
            num_solutions: Number of solutions to generate
            
        Returns:
            List of optimized sequences
        """
        num_codons = len(protein_sequence)
        
        # Create custom problem for expression vs structure
        class ExpressionStructureProblem(CodonOptimizationProblem):
            def _evaluate(self, x, out, *args, **kwargs):
                n_pop = x.shape[0]
                objectives = np.zeros((n_pop, 2))
                
                for i in range(n_pop):
                    dna_sequence = self._decode_solution(x[i])
                    
                    try:
                        optimized_protein = str(Seq(dna_sequence).translate())
                        if optimized_protein != self.original_protein:
                            objectives[i, 0] = 1.0
                            objectives[i, 1] = 1.0
                            continue
                    except:
                        objectives[i, 0] = 1.0
                        objectives[i, 1] = 1.0
                        continue
                    
                    # Objective 1: Negative CAI (maximize expression)
                    cai = self.codon_analyzer.calculate_cai(dna_sequence)
                    objectives[i, 0] = 1.0 - cai
                    
                    # Objective 2: Structure penalty (minimize stable structures)
                    structure_score = self.scorer.calculate_structure_score(dna_sequence)
                    objectives[i, 1] = 1.0 - structure_score
                
                out["F"] = objectives
        
        problem = ExpressionStructureProblem(
            original_protein=protein_sequence,
            scorer=self.scorer,
            codon_analyzer=self.codon_analyzer,
            num_codons=num_codons
        )
        
        algorithm = NSGA2(pop_size=num_solutions * 2)
        res = minimize(problem, algorithm, ('n_gen', 50), verbose=False, seed=42)
        
        solutions = []
        for i, x in enumerate(res.X):
            dna_sequence = problem._decode_solution(x)
            
            try:
                optimized_protein = str(Seq(dna_sequence).translate())
                if optimized_protein != original_protein:
                    continue
            except:
                continue
            
            scores = self.scorer.calculate_total_score(dna_sequence, original_protein)
            
            solution = OptimizedSequence(
                dna_sequence=dna_sequence,
                original_sequence=None,  # No original DNA, only protein
                sequence_id=f"expr_struct_{i+1}",
                metadata={
                    'fitness_score': scores['total_score'],
                    'cai': self.codon_analyzer.calculate_cai(dna_sequence),
                    'structure_score': scores['structure_score'],
                    'pareto_rank': i + 1
                }
            )
            solutions.append(solution)
        
        solutions.sort(key=lambda s: s.metadata.get('fitness_score', 0.0), reverse=True)
        return solutions[:num_solutions]



