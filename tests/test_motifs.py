"""
Tests for motif detection.
"""

import unittest
from codon_optimizer.analysis.motifs import MotifDetector


class TestMotifDetector(unittest.TestCase):
    """Test cases for MotifDetector."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.detector = MotifDetector()
    
    def test_find_tata_boxes(self):
        """Test TATA box detection."""
        sequence = "TATAAA" + "ATCG" * 100
        tata_boxes = self.detector.find_tata_boxes(sequence)
        self.assertGreater(len(tata_boxes), 0)
    
    def test_find_homopolymers(self):
        """Test homopolymer detection."""
        sequence = "AAAAA" + "ATCG" * 100  # Homopolymer A
        homopolymers = self.detector.find_homopolymers(sequence, min_length=5)
        self.assertGreater(len(homopolymers), 0)
    
    def test_find_tandem_repeats(self):
        """Test tandem repeat detection."""
        sequence = "ATCGATCGATCG" + "GCTA" * 50  # Tandem repeat
        repeats = self.detector.find_tandem_repeats(sequence)
        self.assertIsInstance(repeats, list)
    
    def test_calculate_motif_penalty(self):
        """Test motif penalty calculation."""
        sequence = "ATCG" * 100  # Clean sequence
        penalty = self.detector.calculate_motif_penalty(sequence)
        self.assertGreaterEqual(penalty, 0.0)
        self.assertLessEqual(penalty, 1.0)


if __name__ == '__main__':
    unittest.main()




