"""
Antibody-specific optimization handler for heavy and light chains.
"""

from typing import List, Tuple, Dict, Optional
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.core.optimizer import CodonOptimizer
from codon_optimizer.core.pareto import ParetoOptimizer
from codon_optimizer.core.scorer import MultiCriteriaScorer
from codon_optimizer.core.sequence import OptimizedSequence, SequencePair
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from Bio.Seq import Seq


class AntibodyOptimizer:
    """Optimizer for antibody heavy and light chains."""
    
    def __init__(self, config: OptimizationConfig):
        """
        Initialize antibody optimizer.
        
        Args:
            config: Optimization configuration
        """
        self.config = config
        self.scorer = MultiCriteriaScorer(config)
        self.codon_analyzer = CodonUsageAnalyzer(
            host=config.host,
            max_codon_usage=config.max_codon_usage
        )
    
    def optimize_heavy_chain(self, heavy_chain_protein: str) -> OptimizedSequence:
        """
        Optimize heavy chain sequence from protein.
        
        Args:
            heavy_chain_protein: Heavy chain protein sequence (amino acids)
            
        Returns:
            Optimized heavy chain sequence
        """
        optimizer = CodonOptimizer(self.config)
        return optimizer.optimize(heavy_chain_protein)
    
    def optimize_light_chain(self, light_chain_protein: str) -> OptimizedSequence:
        """
        Optimize light chain sequence from protein.
        
        Args:
            light_chain_protein: Light chain protein sequence (amino acids)
            
        Returns:
            Optimized light chain sequence
        """
        optimizer = CodonOptimizer(self.config)
        return optimizer.optimize(light_chain_protein)
    
    def optimize_coordinated(self,
                            heavy_chain_protein: str,
                            light_chain_protein: str) -> SequencePair:
        """
        Optimize heavy and light chains with coordination.
        
        Args:
            heavy_chain_protein: Heavy chain protein sequence (amino acids)
            light_chain_protein: Light chain protein sequence (amino acids)
            
        Returns:
            SequencePair with optimized chains
        """
        # Optimize chains separately first
        hc_optimized = self.optimize_heavy_chain(heavy_chain_protein)
        lc_optimized = self.optimize_light_chain(light_chain_protein)
        
        # Balance expression levels
        if self.config.balance_chains:
            hc_optimized, lc_optimized = self._balance_expression(
                hc_optimized, lc_optimized
            )
        
        return SequencePair(hc_optimized, lc_optimized, "HC_LC")
    
    def _balance_expression(self,
                          hc_sequence: OptimizedSequence,
                          lc_sequence: OptimizedSequence) -> Tuple[OptimizedSequence, OptimizedSequence]:
        """
        Balance expression levels between heavy and light chains.
        
        Args:
            hc_sequence: Heavy chain sequence
            lc_sequence: Light chain sequence
            
        Returns:
            Tuple of (balanced_hc, balanced_lc)
        """
        # Calculate CAI for each chain
        hc_cai = self.codon_analyzer.calculate_cai(hc_sequence.dna_sequence)
        lc_cai = self.codon_analyzer.calculate_cai(lc_sequence.dna_sequence)
        
        # Target: similar CAI values
        target_cai = (hc_cai + lc_cai) / 2
        
        # Adjust if difference is too large
        cai_diff = abs(hc_cai - lc_cai)
        if cai_diff > 0.1:  # Threshold for balancing
            # Adjust the chain with higher CAI to reduce it slightly
            if hc_cai > lc_cai:
                # Reduce HC CAI by using slightly less optimal codons
                hc_sequence = self._adjust_cai(hc_sequence, target_cai, direction='down')
            else:
                # Reduce LC CAI
                lc_sequence = self._adjust_cai(lc_sequence, target_cai, direction='down')
        
        return hc_sequence, lc_sequence
    
    def _adjust_cai(self,
                   sequence: OptimizedSequence,
                   target_cai: float,
                   direction: str = 'down') -> OptimizedSequence:
        """
        Adjust CAI of sequence by modifying codons.
        
        Args:
            sequence: Sequence to adjust
            target_cai: Target CAI value
            direction: 'up' or 'down'
            
        Returns:
            Adjusted sequence
        """
        current_cai = self.codon_analyzer.calculate_cai(sequence.dna_sequence)
        protein = sequence.get_protein_sequence()
        
        if direction == 'down' and current_cai <= target_cai:
            return sequence
        if direction == 'up' and current_cai >= target_cai:
            return sequence
        
        # Modify codons to adjust CAI
        codons = sequence.get_codons()
        modified_codons = []
        
        for i, (codon, aa) in enumerate(zip(codons, protein)):
            if direction == 'down' and current_cai > target_cai:
                # Use less optimal codons
                options = self.codon_analyzer.get_codon_options(aa, avoid_overuse=True)
                if len(options) > 1:
                    # Prefer less optimal options
                    options_sorted = sorted(
                        options,
                        key=lambda c: self.codon_analyzer.codon_table.get(c, 0.0),
                        reverse=True
                    )
                    # Use a middle option
                    new_codon = options_sorted[min(1, len(options_sorted) - 1)]
                    modified_codons.append(new_codon)
                else:
                    modified_codons.append(codon)
            else:
                modified_codons.append(codon)
        
        new_dna = ''.join(modified_codons)
        return OptimizedSequence(
            dna_sequence=new_dna,
            original_sequence=sequence.original_sequence,
            sequence_id=sequence.sequence_id,
            metadata=sequence.metadata
        )
    
    def optimize_bispecific(self,
                           heavy_chain_1_dna: str,
                           light_chain_1_dna: str,
                           heavy_chain_2_dna: str,
                           light_chain_2_dna: str) -> Dict[str, SequencePair]:
        """
        Optimize bispecific antibody with two heavy-light chain pairs.
        
        Args:
            heavy_chain_1_dna: First heavy chain DNA
            light_chain_1_dna: First light chain DNA
            heavy_chain_2_dna: Second heavy chain DNA
            light_chain_2_dna: Second light chain DNA
            
        Returns:
            Dictionary with optimized pairs
        """
        pair1 = self.optimize_coordinated(heavy_chain_1_dna, light_chain_1_dna)
        pair2 = self.optimize_coordinated(heavy_chain_2_dna, light_chain_2_dna)
        
        # Balance between pairs
        hc1_cai = self.codon_analyzer.calculate_cai(pair1.sequence1.dna_sequence)
        hc2_cai = self.codon_analyzer.calculate_cai(pair2.sequence1.dna_sequence)
        
        # Ensure similar expression levels
        if abs(hc1_cai - hc2_cai) > 0.1:
            target = (hc1_cai + hc2_cai) / 2
            if hc1_cai > hc2_cai:
                pair1.sequence1, _ = self._balance_expression(
                    pair1.sequence1, pair1.sequence2
                )
            else:
                pair2.sequence1, _ = self._balance_expression(
                    pair2.sequence1, pair2.sequence2
                )
        
        return {
            'pair1': pair1,
            'pair2': pair2
        }
    
    def transfer_variable_domains(self,
                                  variable_domain_dna: str,
                                  constant_domain_dna: str) -> OptimizedSequence:
        """
        Transfer variable domains to constant domains with codon optimization.
        
        Args:
            variable_domain_dna: Variable domain DNA sequence
            constant_domain_dna: Constant domain DNA sequence
            
        Returns:
            Optimized combined sequence
        """
        # Combine sequences
        combined_dna = variable_domain_dna + constant_domain_dna
        
        # Optimize the combined sequence
        optimizer = CodonOptimizer(self.config)
        optimized = optimizer.optimize(combined_dna)
        
        return optimized
    
    def check_signal_peptide_compatibility(self,
                                         hc_sequence: str,
                                         lc_sequence: str) -> Dict[str, any]:
        """
        Check compatibility of signal peptides between chains.
        
        Args:
            hc_sequence: Heavy chain DNA sequence
            lc_sequence: Light chain DNA sequence
            
        Returns:
            Dictionary with compatibility analysis
        """
        # Extract signal peptides (first ~60 nucleotides typically)
        hc_signal = hc_sequence[:60] if len(hc_sequence) >= 60 else hc_sequence
        lc_signal = lc_sequence[:60] if len(lc_sequence) >= 60 else lc_sequence
        
        # Check for similar structure/composition
        from codon_optimizer.analysis.gc_content import GCContentAnalyzer
        gc_analyzer = GCContentAnalyzer()
        
        hc_gc = gc_analyzer.calculate_gc_content(hc_signal)
        lc_gc = gc_analyzer.calculate_gc_content(lc_signal)
        
        gc_diff = abs(hc_gc - lc_gc)
        
        return {
            'hc_signal_gc': hc_gc,
            'lc_signal_gc': lc_gc,
            'gc_difference': gc_diff,
            'is_compatible': gc_diff < 0.2,  # Threshold for compatibility
            'hc_signal_length': len(hc_signal),
            'lc_signal_length': len(lc_signal)
        }
    
    def get_chain_statistics(self, pair: SequencePair) -> Dict[str, any]:
        """
        Get statistics for heavy and light chain pair.
        
        Args:
            pair: SequencePair object
            
        Returns:
            Dictionary with chain statistics
        """
        hc_cai = self.codon_analyzer.calculate_cai(pair.sequence1.dna_sequence)
        lc_cai = self.codon_analyzer.calculate_cai(pair.sequence2.dna_sequence)
        
        hc_scores = self.scorer.calculate_total_score(
            pair.sequence1.dna_sequence,
            pair.sequence1.get_protein_sequence()
        )
        lc_scores = self.scorer.calculate_total_score(
            pair.sequence2.dna_sequence,
            pair.sequence2.get_protein_sequence()
        )
        
        return {
            'heavy_chain': {
                'cai': hc_cai,
                'length': len(pair.sequence1),
                'scores': hc_scores
            },
            'light_chain': {
                'cai': lc_cai,
                'length': len(pair.sequence2),
                'scores': lc_scores
            },
            'cai_difference': abs(hc_cai - lc_cai),
            'is_balanced': abs(hc_cai - lc_cai) < 0.1
        }



