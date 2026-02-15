"""
Tests for GC content analysis.
"""

import unittest
from codon_optimizer.analysis.gc_content import GCContentAnalyzer


class TestGCContentAnalyzer(unittest.TestCase):
    """Test cases for GCContentAnalyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = GCContentAnalyzer()
    
    def test_calculate_gc_content(self):
        """Test GC content calculation."""
        sequence = "ATCGATCG"  # 50% GC
        gc = self.analyzer.calculate_gc_content(sequence)
        self.assertEqual(gc, 0.5)
    
    def test_calculate_gc_profile(self):
        """Test GC profile calculation."""
        sequence = "ATCGATCGATCG" * 10  # Repeat pattern
        profile = self.analyzer.calculate_gc_profile(sequence)
        self.assertGreater(len(profile), 0)
        self.assertTrue(all(0.0 <= gc <= 1.0 for gc in profile))
    
    def test_find_gc_violations(self):
        """Test finding GC violations."""
        # Create sequence with extreme GC
        sequence = "AAAAA" * 20  # Very low GC
        violations = self.analyzer.find_gc_violations(sequence)
        self.assertIsInstance(violations, list)
    
    def test_calculate_gc_score(self):
        """Test GC score calculation."""
        sequence = "ATCGATCG" * 50  # Balanced GC
        score = self.analyzer.calculate_gc_score(sequence)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)


if __name__ == '__main__':
    unittest.main()




