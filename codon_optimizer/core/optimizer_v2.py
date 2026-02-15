"""
New CodonOptimizer using DnaChisel-inspired architecture.
Maintains backward compatibility with simple interface.
"""

from typing import Optional, List
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.core.dna_optimization_problem import DnaOptimizationProblem, NoSolutionError
from codon_optimizer.core.sequence import OptimizedSequence
from codon_optimizer.specifications import (
    EnforceTranslation,
    EnforceGCContent,
    AvoidPattern,
    AvoidMotifs,
    EnforceRestrictionSites,
    EnforceProteinIntegrity,
    MaximizeCAI,
    OptimizeCodonUsage,
    OptimizeGCContent,
    MinimizeMotifs,
    OptimizeCodonPairs,
    OptimizeMRNAStructure,
)
from codon_optimizer.specifications.specification import Location


class CodonOptimizerV2:
    """
    New codon optimizer with constraint/objective separation.
    
    IMPORTANT: This optimizer works with AMINO ACID SEQUENCES as input.
    All DNA sequences are generated from the protein sequence, ensuring
    that all variants encode the same amino acid sequence (synonymous mutations only).
    
    Maintains simple interface while using DnaChisel-inspired architecture internally.
    """
    
    def __init__(self, config: OptimizationConfig):
        """
        Initialize optimizer.
        
        Args:
            config: Optimization configuration
        """
        self.config = config
    
    def _create_constraints(self, protein_sequence: str) -> List:
        """Create constraint specifications from config."""
        constraints = []
        
        # Enforce protein integrity
        constraints.append(EnforceProteinIntegrity(
            original_protein=protein_sequence,
            priority=100  # Highest priority
        ))
        
        # Enforce translation
        constraints.append(EnforceTranslation(
            protein_sequence=protein_sequence,
            priority=90
        ))
        
        # GC content constraint
        constraints.append(EnforceGCContent(
            min_gc=self.config.gc_min,
            max_gc=self.config.gc_max,
            window=self.config.gc_window_size,
            priority=50
        ))
        
        # Restriction sites
        if self.config.restriction_sites_to_remove:
            constraints.append(EnforceRestrictionSites(
                sites_to_remove=self.config.restriction_sites_to_remove,
                sites_to_keep=self.config.restriction_sites_to_keep,
                priority=40
            ))
        
        # Avoid problematic motifs (soft constraint - will be optimized in objectives)
        # For long sequences, make it an objective rather than hard constraint
        if len(protein_sequence) > 100:
            # For long sequences, minimize motifs as objective, not hard constraint
            pass  # Will be handled in objectives
        else:
            constraints.append(AvoidMotifs(priority=30))
        
        return constraints
    
    def _create_objectives(self) -> List:
        """Create objective specifications from config."""
        objectives = []
        
        # Maximize CAI
        objectives.append(MaximizeCAI(
            host=self.config.host,
            target_min=self.config.target_cai_min,
            target_max=self.config.target_cai_max,
            boost=self.config.weight_cai
        ))
        
        # Optimize codon usage (avoid tRNA depletion)
        objectives.append(OptimizeCodonUsage(
            host=self.config.host,
            max_codon_usage=self.config.max_codon_usage,
            boost=self.config.weight_codon_pairs
        ))
        
        # Optimize GC content
        target_gc = (self.config.gc_min + self.config.gc_max) / 2
        objectives.append(OptimizeGCContent(
            target_gc=target_gc,
            tolerance=(self.config.gc_max - self.config.gc_min) / 2,
            boost=self.config.weight_gc
        ))
        
        # Minimize motifs
        objectives.append(MinimizeMotifs(boost=self.config.weight_motifs))
        
        # Optimize codon pairs
        objectives.append(OptimizeCodonPairs(
            host=self.config.host,
            boost=self.config.weight_codon_pairs
        ))
        
        # Optimize mRNA structure
        objectives.append(OptimizeMRNAStructure(
            max_5prime_dg=self.config.max_5prime_dg,
            boost=self.config.weight_mrna_structure
        ))
        
        return objectives
    
    def optimize(self,
                protein_sequence: str,
                max_combinations: Optional[int] = None,
                progress_callback=None) -> OptimizedSequence:
        """
        Optimize codon usage from amino acid sequence.
        
        IMPORTANT: Input is an AMINO ACID SEQUENCE (protein).
        All DNA sequences are generated from this protein sequence.
        All mutations are SYNONYMOUS - they preserve the amino acid sequence.
        
        Uses new architecture:
        1. Create problem with constraints and objectives
        2. Generate initial DNA sequence from protein using optimal codons
        3. resolve_constraints() - solve all constraints (using synonymous codon mutations)
        4. optimize() - optimize objectives while maintaining protein sequence
        
        Args:
            protein_sequence: Protein sequence (amino acids, single-letter code)
                             Example: "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
            max_combinations: Ignored (kept for compatibility)
            progress_callback: Optional callback function(message) for progress updates
            
        Returns:
            OptimizedSequence with optimized DNA sequence encoding the same protein
        """
        def log(msg):
            if progress_callback:
                progress_callback(msg)
            else:
                import sys
                print(msg, file=sys.stderr, flush=True)
        
        import time
        overall_start = time.time()
        log(f"[1/5] Creating constraints and objectives...")
        # Create constraints and objectives
        constraints = self._create_constraints(protein_sequence)
        objectives = self._create_objectives()
        log(f"  ✓ Created {len(constraints)} constraints, {len(objectives)} objectives")
        
        log(f"[2/5] Generating initial DNA sequence from protein ({len(protein_sequence)} aa)...")
        # Create optimization problem
        problem = DnaOptimizationProblem(
            sequence="",  # Will be generated from protein
            protein_sequence=protein_sequence,
            constraints=constraints,
            objectives=objectives,
            host=self.config.host
        )
        log(f"  ✓ Initial sequence generated ({len(problem.sequence)} bp)")
        
        # Adjust algorithm parameters based on sequence length
        if len(protein_sequence) > 100:
            problem.randomization_threshold = 5000
            problem.max_random_iters = 2000
            log(f"  ✓ Using parameters for long sequence (threshold: {problem.randomization_threshold})")
        else:
            problem.randomization_threshold = 10000
            problem.max_random_iters = 1000
            log(f"  ✓ Using parameters for short sequence (threshold: {problem.randomization_threshold})")
        
        # Solve constraints
        log(f"[3/5] Resolving constraints...")
        try:
            problem.resolve_constraints(progress_callback=log)
            log(f"  ✓ All constraints resolved successfully")
        except NoSolutionError as e:
            log(f"  ✗ Constraint resolution failed: {e}")
            # Fallback to old method if new method fails
            log(f"[FALLBACK] Using legacy optimizer...")
            from codon_optimizer.core.optimizer import CodonOptimizer
            old_optimizer = CodonOptimizer(self.config)
            return old_optimizer.optimize(protein_sequence, max_combinations)
        
        # Optimize objectives
        log(f"[4/5] Optimizing objectives...")
        try:
            problem.optimize(progress_callback=log)
            log(f"  ✓ Objectives optimization completed")
        except Exception as e:
            log(f"  ⚠ Optimization warning: {e} (continuing with current solution)")
            # If optimization fails, at least we have a valid sequence
            pass
        
        # Get final score
        log(f"[5/5] Calculating final scores...")
        from codon_optimizer.core.scorer import MultiCriteriaScorer
        scorer = MultiCriteriaScorer(self.config)
        score_dict = scorer.calculate_total_score(problem.sequence, protein_sequence)
        total_time = time.time() - overall_start
        log(f"  ✓ Final fitness score: {score_dict['total_score']:.4f}")
        log(f"  ✓ Total optimization time: {total_time:.1f}s")
        
        return OptimizedSequence(
            dna_sequence=problem.sequence,
            original_sequence=None,
            sequence_id="optimized",
            metadata={
                'fitness_score': score_dict['total_score'],
                'optimization_method': 'dna_chisel_inspired',
                'constraints_passed': problem.all_constraints_pass(),
                'objectives_score': problem.objectives_evaluations().scores_sum(),
                **score_dict  # Include all individual scores
            }
        )
    
    def optimize_multiple_diverse(self,
                                  protein_sequence: str,
                                  num_solutions: int = 5,
                                  progress_callback=None) -> List[OptimizedSequence]:
        """
        Generate multiple optimized solutions, each optimized for different criteria.
        
        Generates solutions optimized for:
        1. Best overall score (balanced)
        2. Best CAI (codon adaptation)
        3. Best GC content
        4. Best motif score (minimal problematic motifs)
        5. Best codon pair usage
        
        Args:
            protein_sequence: Protein sequence (amino acids)
            num_solutions: Number of diverse solutions to generate (default: 5)
            progress_callback: Optional callback function(message) for progress updates
            
        Returns:
            List of OptimizedSequence objects, each optimized for different criteria
        """
        def log(msg):
            if progress_callback:
                progress_callback(msg)
            else:
                import sys
                print(msg, file=sys.stderr, flush=True)
        
        import time
        overall_start = time.time()
        
        # Strategy: optimize with different objective weights
        strategies = [
            {
                'name': 'balanced',
                'description': 'Best overall score (balanced)',
                'weights': {
                    'cai': 1.0,
                    'codon_usage': 1.0,
                    'gc': 1.0,
                    'motifs': 1.0,
                    'codon_pairs': 1.0,
                    'mrna_structure': 1.0
                }
            },
            {
                'name': 'high_cai',
                'description': 'Best CAI (codon adaptation)',
                'weights': {
                    'cai': 3.0,
                    'codon_usage': 1.0,
                    'gc': 0.5,
                    'motifs': 0.5,
                    'codon_pairs': 0.5,
                    'mrna_structure': 0.5
                }
            },
            {
                'name': 'optimal_gc',
                'description': 'Best GC content',
                'weights': {
                    'cai': 0.5,
                    'codon_usage': 0.5,
                    'gc': 3.0,
                    'motifs': 0.5,
                    'codon_pairs': 0.5,
                    'mrna_structure': 0.5
                }
            },
            {
                'name': 'minimal_motifs',
                'description': 'Minimal problematic motifs',
                'weights': {
                    'cai': 0.5,
                    'codon_usage': 0.5,
                    'gc': 0.5,
                    'motifs': 3.0,
                    'codon_pairs': 0.5,
                    'mrna_structure': 0.5
                }
            },
            {
                'name': 'best_codon_pairs',
                'description': 'Best codon pair usage',
                'weights': {
                    'cai': 0.5,
                    'codon_usage': 1.0,
                    'gc': 0.5,
                    'motifs': 0.5,
                    'codon_pairs': 3.0,
                    'mrna_structure': 0.5
                }
            }
        ]
        
        # Limit to requested number
        strategies = strategies[:num_solutions]
        
        solutions = []
        from codon_optimizer.core.scorer import MultiCriteriaScorer
        
        for i, strategy in enumerate(strategies, 1):
            log(f"[Strategy {i}/{len(strategies)}] Optimizing for: {strategy['description']}...")
            
            # Create config with modified weights
            strategy_config = OptimizationConfig(
                host=self.config.host,
                weight_cai=strategy['weights']['cai'] * self.config.weight_cai,
                weight_codon_pairs=strategy['weights']['codon_pairs'] * self.config.weight_codon_pairs,
                weight_gc=strategy['weights']['gc'] * self.config.weight_gc,
                weight_motifs=strategy['weights']['motifs'] * self.config.weight_motifs,
                weight_mrna_structure=strategy['weights']['mrna_structure'] * self.config.weight_mrna_structure,
                gc_min=self.config.gc_min,
                gc_max=self.config.gc_max,
                gc_window_size=self.config.gc_window_size,
                target_cai_min=self.config.target_cai_min,
                target_cai_max=self.config.target_cai_max,
                max_codon_usage=self.config.max_codon_usage,
                max_5prime_dg=self.config.max_5prime_dg,
                restriction_sites_to_remove=self.config.restriction_sites_to_remove,
                restriction_sites_to_keep=self.config.restriction_sites_to_keep
            )
            
            # Create optimizer with strategy config
            strategy_optimizer = CodonOptimizerV2(strategy_config)
            
            # Optimize
            optimized = strategy_optimizer.optimize(
                protein_sequence,
                progress_callback=lambda msg: log(f"  {msg}") if progress_callback else None
            )
            
            # Calculate detailed scores
            scorer = MultiCriteriaScorer(strategy_config)
            score_dict = scorer.calculate_total_score(optimized.dna_sequence, protein_sequence)
            
            # Update metadata
            optimized.sequence_id = f"optimized_{strategy['name']}"
            optimized.metadata.update({
                'strategy': strategy['name'],
                'strategy_description': strategy['description'],
                'fitness_score': score_dict['total_score'],
                **score_dict
            })
            
            solutions.append(optimized)
            log(f"  ✓ {strategy['description']}: score={score_dict['total_score']:.4f}, CAI={score_dict.get('cai', 0):.4f}, GC={score_dict.get('gc_score', 0):.4f}")
        
        total_time = time.time() - overall_start
        log(f"✓ Generated {len(solutions)} diverse solutions in {total_time:.1f}s")
        
        return solutions

