"""
Tests for codon usage analysis.
"""

import unittest
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer


class TestCodonUsageAnalyzer(unittest.TestCase):
    """Test cases for CodonUsageAnalyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = CodonUsageAnalyzer(host="CHO")
    
    def test_calculate_cai(self):
        """Test CAI calculation."""
        # Test sequence with optimal codons
        sequence = "ATGGCCACC"  # Met-Ala-Thr (all optimal)
        cai = self.analyzer.calculate_cai(sequence)
        self.assertGreater(cai, 0.0)
        self.assertLessEqual(cai, 1.0)
    
    def test_get_optimal_codon(self):
        """Test getting optimal codon for amino acid."""
        optimal = self.analyzer.get_optimal_codon('M')  # Methionine
        self.assertEqual(optimal, 'ATG')
    
    def test_check_tRNA_depletion(self):
        """Test tRNA depletion check."""
        # Create sequence with repeated codon
        sequence = "ATGGCCGCCGCCGCCGCC"  # Many GCC (Ala)
        has_risk, ratios = self.analyzer.check_tRNA_depletion(sequence)
        # Should detect overuse
        self.assertIsInstance(has_risk, bool)
        self.assertIsInstance(ratios, dict)
    
    def test_get_codon_options(self):
        """Test getting codon options for amino acid."""
        options = self.analyzer.get_codon_options('A')  # Alanine
        self.assertGreater(len(options), 0)
        self.assertIn('GCC', options)  # Should include common codons


if __name__ == '__main__':
    unittest.main()




