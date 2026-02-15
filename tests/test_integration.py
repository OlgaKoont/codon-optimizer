"""
Integration tests for full optimization workflow.
"""

import unittest
import tempfile
from pathlib import Path
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.core.optimizer import CodonOptimizer
from codon_optimizer.io.fasta_handler import FastaHandler
from codon_optimizer.core.scorer import MultiCriteriaScorer


class TestIntegration(unittest.TestCase):
    """Integration tests for full workflow."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.config = OptimizationConfig(
            host="CHO",
            population_size=20,
            max_generations=10
        )
        
        # Create a simple test sequence
        self.test_sequence = "ATGGCCACCATGGCCACC"  # Met-Ala-Thr-Met-Ala-Thr
    
    def test_full_optimization(self):
        """Test full optimization workflow."""
        optimizer = CodonOptimizer(self.config)
        optimized = optimizer.optimize(self.test_sequence)
        
        # Check that protein sequence is preserved
        original_protein = FastaHandler.translate_dna_to_protein(self.test_sequence)
        optimized_protein = optimized.get_protein_sequence()
        self.assertEqual(original_protein, optimized_protein)
        
        # Check that CAI improved or stayed reasonable
        scorer = MultiCriteriaScorer(self.config)
        scores = scorer.calculate_total_score(
            optimized.dna_sequence,
            original_protein
        )
        self.assertGreater(scores['total_score'], 0.0)
    
    def test_fasta_io(self):
        """Test FASTA file I/O."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">test_sequence\n")
            f.write(self.test_sequence + "\n")
            temp_file = f.name
        
        try:
            # Read
            seq_id, sequence = FastaHandler.read_single_sequence(temp_file)
            self.assertEqual(seq_id, "test_sequence")
            self.assertEqual(sequence, self.test_sequence)
            
            # Write
            output_file = temp_file.replace('.fasta', '_out.fasta')
            FastaHandler.write_single_sequence("output", sequence, output_file)
            
            # Read back
            seq_id2, sequence2 = FastaHandler.read_single_sequence(output_file)
            self.assertEqual(sequence2, sequence)
            
            Path(output_file).unlink()
        finally:
            Path(temp_file).unlink()
    
    def test_scorer_integration(self):
        """Test scorer with all analyzers."""
        scorer = MultiCriteriaScorer(self.config)
        original_protein = FastaHandler.translate_dna_to_protein(self.test_sequence)
        
        analysis = scorer.get_detailed_analysis(self.test_sequence, original_protein)
        
        # Check that all components are present
        self.assertIn('scores', analysis)
        self.assertIn('cai', analysis)
        self.assertIn('gc_statistics', analysis)
        self.assertIn('motifs', analysis)
        self.assertIn('restriction_sites', analysis)


if __name__ == '__main__':
    unittest.main()




