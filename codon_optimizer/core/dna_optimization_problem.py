"""
DNA Optimization Problem with constraint/objective separation.
Inspired by DnaChisel architecture.
"""

from typing import List, Optional, Dict
from codon_optimizer.specifications import (
    Constraint, Objective,
    SpecEvaluation, SpecEvaluations
)
from codon_optimizer.core.mutation_space import MutationSpace
from codon_optimizer.core.codon_generator import CodonCombinationGenerator


class NoSolutionError(Exception):
    """Raised when no solution can be found."""
    pass


class DnaOptimizationProblem:
    """
    DNA optimization problem with constraints and objectives.
    
    IMPORTANT: This class works with AMINO ACID SEQUENCES as input.
    All DNA sequences are generated from the protein sequence, ensuring
    that all variants encode the same amino acid sequence (synonymous mutations only).
    
    Architecture:
    1. Define problem with protein_sequence (amino acids), constraints, and objectives
    2. Generate initial DNA sequence from protein using optimal codons
    3. resolve_constraints() - solve all hard constraints (using synonymous codon mutations)
    4. optimize() - optimize objectives while maintaining constraints and protein sequence
    """
    
    # Algorithm parameters
    randomization_threshold = 10000  # Use exhaustive search if space < threshold
    max_random_iters = 1000
    mutations_per_iteration = 2
    optimization_stagnation_tolerance = 100
    local_extensions = (0, 5)  # Extend local zone by these amounts if needed
    
    def __init__(self,
                 sequence: str,
                 protein_sequence: Optional[str] = None,
                 constraints: Optional[List[Constraint]] = None,
                 objectives: Optional[List[Objective]] = None,
                 host: str = "CHO"):
        """
        Initialize optimization problem from amino acid sequence.
        
        IMPORTANT: Provide protein_sequence (amino acids) as input.
        The DNA sequence will be generated automatically, and all mutations
        will preserve the amino acid sequence (synonymous mutations only).
        
        Args:
            sequence: Initial DNA sequence (empty string "" to auto-generate from protein)
            protein_sequence: Protein sequence (amino acids, single-letter code) - REQUIRED
            constraints: List of constraint specifications
            objectives: List of objective specifications
            host: Host organism (for codon usage preferences)
            
        Raises:
            ValueError: If neither sequence nor protein_sequence provided,
                      or if sequence doesn't translate to protein_sequence
        """
        self.host = host
        self.codon_generator = CodonCombinationGenerator(host=host)
        
        # Generate sequence from protein if needed
        if protein_sequence and not sequence:
            sequence = self._generate_initial_sequence(protein_sequence)
        elif protein_sequence and sequence:
            # Validate that sequence translates to protein
            from Bio.Seq import Seq
            if str(Seq(sequence).translate()) != protein_sequence:
                raise ValueError("Sequence does not translate to provided protein")
        
        if not sequence:
            raise ValueError("Either sequence or protein_sequence must be provided")
        
        self.sequence = sequence.upper()
        self.sequence_before = self.sequence
        self.protein_sequence = protein_sequence
        
        # Initialize constraints and objectives
        self.constraints = constraints or []
        self.objectives = objectives or []
        
        # Create mutation space if protein sequence provided
        if protein_sequence:
            self.mutation_space = MutationSpace(protein_sequence, host=host)
        else:
            # For DNA-only sequences, mutation space is limited
            self.mutation_space = None
        
        # Evaluation cache
        self._constraints_before = None
        self._objectives_before = None
    
    def _generate_initial_sequence(self, protein_sequence: str) -> str:
        """Generate initial sequence using optimal codons."""
        from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
        analyzer = CodonUsageAnalyzer(host=self.host)
        return ''.join(analyzer.get_optimal_codon(aa) for aa in protein_sequence)
    
    def constraints_evaluations(self) -> SpecEvaluations:
        """Evaluate all constraints."""
        evaluations = []
        for constraint in self.constraints:
            eval_ = constraint.evaluate(self.sequence, self.protein_sequence)
            evaluations.append((constraint, eval_))
        return SpecEvaluations(evaluations)
    
    def objectives_evaluations(self) -> SpecEvaluations:
        """Evaluate all objectives."""
        evaluations = []
        for objective in self.objectives:
            eval_ = objective.evaluate(self.sequence, self.protein_sequence)
            evaluations.append((objective, eval_))
        return SpecEvaluations(evaluations)
    
    def all_constraints_pass(self) -> bool:
        """Check if all constraints pass."""
        if not self.constraints:
            return True
        evaluations = self.constraints_evaluations()
        return all(eval_.passes for _, eval_ in evaluations.evaluations)
    
    def constraints_text_summary(self, failed_only: bool = False) -> str:
        """Get text summary of constraints."""
        evaluations = self.constraints_evaluations()
        if failed_only:
            evaluations = evaluations.filter('failing')
        return evaluations.to_text()
    
    def objectives_text_summary(self) -> str:
        """Get text summary of objectives."""
        evaluations = self.objectives_evaluations()
        return evaluations.to_text()
    
    def resolve_constraints(self, progress_callback=None):
        """
        Resolve all constraints.
        
        Uses local iterative approach:
        1. Find constraint violations
        2. For each violation, create local problem
        3. Solve locally
        4. Repeat until all constraints pass
        """
        def log(msg):
            if progress_callback:
                progress_callback(msg)
            else:
                import sys
                print(msg, flush=True)
        
        if not self.constraints:
            log("  No constraints to resolve")
            return
        
        # Sort constraints by priority
        constraints_to_solve = sorted(
            [c for c in self.constraints if not getattr(c, 'enforced_by_nucleotide_restrictions', False)],
            key=lambda c: -c.priority
        )
        
        log(f"  Resolving {len(constraints_to_solve)} constraints (sorted by priority)...")
        import time
        constraints_start = time.time()
        
        for i, constraint in enumerate(constraints_to_solve, 1):
            constraint_name = constraint.__class__.__name__
            elapsed = time.time() - constraints_start
            log(f"  [{i}/{len(constraints_to_solve)}] Checking {constraint_name}... (elapsed: {elapsed:.1f}s)")
            constraint_start = time.time()
            self._resolve_constraint(constraint, progress_callback=log)
            constraint_time = time.time() - constraint_start
            log(f"    ✓ {constraint_name} resolved (took {constraint_time:.1f}s)")
        
        # Final check
        if not self.all_constraints_pass():
            failed = self.constraints_text_summary(failed_only=True)
            raise NoSolutionError(
                f"Failed to resolve all constraints:\n{failed}"
            )
        log(f"  ✓ All constraints satisfied")
    
    def _resolve_constraint(self, constraint: Constraint, progress_callback=None):
        """Resolve a single constraint."""
        def log(msg):
            if progress_callback:
                progress_callback(f"    {msg}")
            else:
                import sys
                print(f"    {msg}", flush=True)
        
        evaluation = constraint.evaluate(self.sequence, self.protein_sequence)
        if evaluation.passes:
            return
        
        log(f"Constraint violated: {evaluation.message}")
        
        # If mutation space is available, use it
        if self.mutation_space and self.mutation_space.space_size < self.randomization_threshold:
            # Exhaustive search
            log(f"Using exhaustive search (space size: {self.mutation_space.space_size})")
            self._resolve_constraint_by_exhaustive_search(constraint, progress_callback=log)
        else:
            # Random mutations
            space_size = self.mutation_space.space_size if self.mutation_space else "unknown"
            log(f"Using random mutations (space size: {space_size}, max iterations: {self.max_random_iters})")
            self._resolve_constraint_by_random_mutations(constraint, progress_callback=log)
    
    def _resolve_constraint_by_exhaustive_search(self, constraint: Constraint, progress_callback=None):
        """Resolve constraint by trying all variants."""
        import time
        
        if not self.mutation_space:
            raise NoSolutionError("Cannot use exhaustive search without mutation space")
        
        total = self.mutation_space.space_size
        checked = 0
        start_time = time.time()
        last_update_time = start_time
        update_interval = 2.0  # Update every 2 seconds
        
        for variant in self.mutation_space.all_variants():
            checked += 1
            current_time = time.time()
            elapsed = current_time - start_time
            
            # Update every 2 seconds OR every 100 variants OR every 10% progress
            should_update = (
                (current_time - last_update_time >= update_interval) or
                (checked % 100 == 0) or
                (checked % max(1, total // 10) == 0)
            )
            
            if should_update and progress_callback:
                rate = checked / elapsed if elapsed > 0 else 0
                eta = (total - checked) / rate if rate > 0 else 0
                progress_callback(
                    f"      Checking variant {checked}/{total} ({100*checked//total:.1f}%, "
                    f"elapsed: {elapsed:.1f}s, rate: {rate:.1f} var/s, ETA: {eta:.1f}s)..."
                )
                last_update_time = current_time
            
            self.sequence = variant
            if constraint.evaluate(self.sequence, self.protein_sequence).passes:
                if self.all_constraints_pass():
                    if progress_callback:
                        elapsed = time.time() - start_time
                        progress_callback(f"      ✓ Solution found at variant {checked}/{total} (took {elapsed:.1f}s)")
                    return
        
        raise NoSolutionError(f"Exhaustive search failed for {constraint}")
    
    def _resolve_constraint_by_random_mutations(self, constraint: Constraint, progress_callback=None):
        """Resolve constraint by random mutations."""
        import time
        
        if not self.mutation_space:
            raise NoSolutionError("Cannot use random mutations without mutation space")
        
        evaluation = constraint.evaluate(self.sequence, self.protein_sequence)
        score = evaluation.score if not evaluation.passes else 1.0
        
        start_time = time.time()
        last_update_time = start_time
        update_interval = 2.0  # Update every 2 seconds
        
        for i in range(self.max_random_iters):
            current_time = time.time()
            elapsed = current_time - start_time
            
            # Update every 2 seconds OR every 50 iterations OR every 10% progress
            should_update = (
                (current_time - last_update_time >= update_interval) or
                (i % 50 == 0) or
                (i % max(1, self.max_random_iters // 10) == 0)
            )
            
            if should_update and progress_callback:
                rate = i / elapsed if elapsed > 0 else 0
                eta = (self.max_random_iters - i) / rate if rate > 0 else 0
                progress_callback(
                    f"      Iteration {i}/{self.max_random_iters} "
                    f"(score: {score:.3f}, elapsed: {elapsed:.1f}s, "
                    f"rate: {rate:.1f} iter/s, ETA: {eta:.1f}s)..."
                )
                last_update_time = current_time
            
            if constraint.evaluate(self.sequence, self.protein_sequence).passes:
                if self.all_constraints_pass():
                    if progress_callback:
                        elapsed = time.time() - start_time
                        progress_callback(f"      ✓ Solution found at iteration {i+1} (took {elapsed:.1f}s)")
                    return
            
            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                self.sequence,
                n_mutations=self.mutations_per_iteration
            )
            
            new_evaluation = constraint.evaluate(self.sequence, self.protein_sequence)
            # Use score directly, ensuring it's never 0.0
            new_score = max(0.1, new_evaluation.score)
            
            # Accept if score improved OR if constraint now passes
            if new_score > score or (new_evaluation.passes and not evaluation.passes):
                score = new_score
                evaluation = new_evaluation
                if progress_callback and new_score > score:
                    progress_callback(f"      ↻ Score improved: {score:.3f} (iteration {i+1})")
                if progress_callback and new_evaluation.passes:
                    progress_callback(f"      ✓ Constraint now passes! (iteration {i+1})")
            else:
                self.sequence = previous_sequence
        
        raise NoSolutionError(
            f"Random search failed for {constraint} after {self.max_random_iters} iterations"
        )
    
    def optimize(self, progress_callback=None):
        """
        Optimize objectives while maintaining constraints.
        
        Only optimizes if all constraints pass.
        """
        def log(msg):
            if progress_callback:
                progress_callback(msg)
            else:
                import sys
                print(msg, flush=True)
        
        if not self.all_constraints_pass():
            summary = self.constraints_text_summary(failed_only=True)
            raise ValueError(
                f"Cannot optimize: constraints not satisfied:\n{summary}"
            )
        
        if not self.objectives:
            log("  No objectives to optimize")
            return
        
        # Use appropriate optimization strategy
        if self.mutation_space and self.mutation_space.space_size < self.randomization_threshold:
            log(f"  Using exhaustive search (space size: {self.mutation_space.space_size})")
            self._optimize_by_exhaustive_search(progress_callback=log)
        else:
            space_size = self.mutation_space.space_size if self.mutation_space else "unknown"
            log(f"  Using random mutations (space size: {space_size}, max iterations: {self.max_random_iters})")
            self._optimize_by_random_mutations(progress_callback=log)
    
    def _optimize_by_exhaustive_search(self, progress_callback=None):
        """Optimize by trying all variants."""
        import time
        
        if not self.mutation_space:
            return
        
        current_best_score = self.objectives_evaluations().scores_sum()
        current_best_sequence = self.sequence
        
        total = self.mutation_space.space_size
        checked = 0
        start_time = time.time()
        last_update_time = start_time
        update_interval = 2.0  # Update every 2 seconds
        
        for variant in self.mutation_space.all_variants():
            checked += 1
            current_time = time.time()
            elapsed = current_time - start_time
            
            # Update every 2 seconds OR every 100 variants OR every 10% progress
            should_update = (
                (current_time - last_update_time >= update_interval) or
                (checked % 100 == 0) or
                (checked % max(1, total // 10) == 0)
            )
            
            if should_update and progress_callback:
                rate = checked / elapsed if elapsed > 0 else 0
                eta = (total - checked) / rate if rate > 0 else 0
                progress_callback(
                    f"    Checking variant {checked}/{total} ({100*checked//total:.1f}%, "
                    f"elapsed: {elapsed:.1f}s, rate: {rate:.1f} var/s, ETA: {eta:.1f}s, "
                    f"best: {current_best_score:.4f})..."
                )
                last_update_time = current_time
            
            self.sequence = variant
            if self.all_constraints_pass():
                score = self.objectives_evaluations().scores_sum()
                if score > current_best_score:
                    current_best_score = score
                    current_best_sequence = self.sequence
                    if progress_callback:
                        elapsed = time.time() - start_time
                        progress_callback(
                            f"    ✓ New best score: {score:.4f} at variant {checked} "
                            f"(elapsed: {elapsed:.1f}s)"
                        )
        
        self.sequence = current_best_sequence
        if progress_callback:
            elapsed = time.time() - start_time
            progress_callback(
                f"  ✓ Optimization complete (best score: {current_best_score:.4f}, "
                f"total time: {elapsed:.1f}s)"
            )
    
    def _optimize_by_random_mutations(self, progress_callback=None):
        """Optimize by random mutations."""
        import time
        
        if not self.mutation_space:
            return
        
        score = self.objectives_evaluations().scores_sum()
        stagnating_iterations = 0
        best_score = score
        
        start_time = time.time()
        last_update_time = start_time
        update_interval = 2.0  # Update every 2 seconds
        
        for iteration in range(self.max_random_iters):
            # Check stagnation
            if (self.optimization_stagnation_tolerance is not None and
                stagnating_iterations > self.optimization_stagnation_tolerance):
                if progress_callback:
                    elapsed = time.time() - start_time
                    progress_callback(f"  Stopping due to stagnation after {iteration} iterations ({elapsed:.1f}s)")
                break
            
            current_time = time.time()
            elapsed = current_time - start_time
            
            # Update every 2 seconds OR every 50 iterations OR every 10% progress
            should_update = (
                (current_time - last_update_time >= update_interval) or
                (iteration % 50 == 0) or
                (iteration % max(1, self.max_random_iters // 10) == 0)
            )
            
            if should_update and progress_callback:
                rate = iteration / elapsed if elapsed > 0 else 0
                eta = (self.max_random_iters - iteration) / rate if rate > 0 else 0
                progress_callback(
                    f"    Iteration {iteration}/{self.max_random_iters} "
                    f"(score: {score:.4f}, best: {best_score:.4f}, "
                    f"elapsed: {elapsed:.1f}s, rate: {rate:.1f} iter/s, ETA: {eta:.1f}s)..."
                )
                last_update_time = current_time
            
            previous_sequence = self.sequence
            self.sequence = self.mutation_space.apply_random_mutations(
                self.sequence,
                n_mutations=self.mutations_per_iteration
            )
            
            if self.all_constraints_pass():
                new_score = self.objectives_evaluations().scores_sum()
                if new_score > score:
                    score = new_score
                    if new_score > best_score:
                        best_score = new_score
                        if progress_callback:
                            elapsed = time.time() - start_time
                            progress_callback(
                                f"    ✓ New best score: {best_score:.4f} at iteration {iteration+1} "
                                f"(elapsed: {elapsed:.1f}s)"
                            )
                    stagnating_iterations = 0
                else:
                    self.sequence = previous_sequence
                    stagnating_iterations += 1
            else:
                self.sequence = previous_sequence
                stagnating_iterations += 1
        
        if progress_callback:
            elapsed = time.time() - start_time
            progress_callback(
                f"  ✓ Optimization complete (final score: {best_score:.4f}, "
                f"iterations: {iteration+1}, total time: {elapsed:.1f}s)"
            )

