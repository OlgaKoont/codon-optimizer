"""
GC content analysis and smoothing for DNA sequences.
"""

from typing import List, Tuple, Dict
import numpy as np
from scipy import signal


class GCContentAnalyzer:
    """Analyzer for GC content with smoothing capabilities."""
    
    def __init__(self, 
                 gc_min: float = 0.30,
                 gc_max: float = 0.80,
                 window_size: int = 50):
        """
        Initialize GC content analyzer.
        
        Args:
            gc_min: Minimum allowed GC content (0-1)
            gc_max: Maximum allowed GC content (0-1)
            window_size: Size of sliding window for GC analysis
        """
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.window_size = window_size
    
    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate overall GC content of sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            GC content as fraction (0-1)
        """
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        total = len(sequence)
        
        if total == 0:
            return 0.0
        
        return gc_count / total
    
    def calculate_gc_profile(self, sequence: str) -> List[float]:
        """
        Calculate GC content profile using sliding window.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            List of GC content values for each window position
        """
        sequence = sequence.upper()
        gc_profile = []
        
        for i in range(len(sequence) - self.window_size + 1):
            window = sequence[i:i + self.window_size]
            gc = self.calculate_gc_content(window)
            gc_profile.append(gc)
        
        return gc_profile
    
    def find_gc_violations(self, sequence: str) -> List[Tuple[int, int, float]]:
        """
        Find regions with GC content outside allowed range.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            List of (start_pos, end_pos, gc_content) tuples for violating regions
        """
        gc_profile = self.calculate_gc_profile(sequence)
        violations = []
        
        for i, gc in enumerate(gc_profile):
            if gc < self.gc_min or gc > self.gc_max:
                start_pos = i
                end_pos = i + self.window_size
                violations.append((start_pos, end_pos, gc))
        
        return violations
    
    def smooth_gc_profile(self, 
                          gc_profile: List[float],
                          method: str = "moving_average",
                          window: int = 5) -> List[float]:
        """
        Smooth GC profile to reduce sharp spikes.
        
        Args:
            gc_profile: List of GC content values
            method: Smoothing method ("moving_average" or "savitzky_golay")
            window: Smoothing window size
            
        Returns:
            Smoothed GC profile
        """
        if len(gc_profile) < window:
            return gc_profile
        
        if method == "moving_average":
            smoothed = np.convolve(gc_profile, np.ones(window)/window, mode='same')
        elif method == "savitzky_golay":
            # Use Savitzky-Golay filter for better preservation of features
            smoothed = signal.savgol_filter(gc_profile, window, 3)
        else:
            raise ValueError(f"Unknown smoothing method: {method}")
        
        return smoothed.tolist()
    
    def calculate_gc_score(self, sequence: str) -> float:
        """
        Calculate GC content score (0-1, higher is better).
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Score between 0 and 1
        """
        gc_profile = self.calculate_gc_profile(sequence)
        
        if not gc_profile:
            return 0.0
        
        # Calculate penalty for each violation
        violations = 0
        total_penalty = 0.0
        
        for gc in gc_profile:
            if gc < self.gc_min:
                # Penalty for too low GC
                penalty = (self.gc_min - gc) / self.gc_min
                total_penalty += penalty
                violations += 1
            elif gc > self.gc_max:
                # Penalty for too high GC
                penalty = (gc - self.gc_max) / (1.0 - self.gc_max)
                total_penalty += penalty
                violations += 1
        
        # Normalize penalty
        if len(gc_profile) > 0:
            normalized_penalty = total_penalty / len(gc_profile)
        else:
            normalized_penalty = 1.0
        
        # Score is inverse of penalty
        score = max(0.0, 1.0 - normalized_penalty)
        
        return score
    
    def calculate_gc_variance(self, sequence: str) -> float:
        """
        Calculate variance in GC content along sequence.
        Lower variance indicates smoother GC profile.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Variance of GC content
        """
        gc_profile = self.calculate_gc_profile(sequence)
        
        if len(gc_profile) < 2:
            return 0.0
        
        return float(np.var(gc_profile))
    
    def suggest_codon_changes_for_gc(self,
                                    sequence: str,
                                    target_gc: float,
                                    position: int) -> List[str]:
        """
        Suggest codon changes to adjust GC content at specific position.
        
        Args:
            sequence: DNA sequence
            target_gc: Target GC content
            position: Position to adjust (codon position, 0-indexed)
            codon_usage_analyzer: CodonUsageAnalyzer instance
            
        Returns:
            List of suggested alternative codons
        """
        # This will be used by optimizer to suggest changes
        # For now, return empty list - will be integrated with optimizer
        return []
    
    def get_gc_statistics(self, sequence: str) -> Dict[str, float]:
        """
        Get comprehensive GC statistics for sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with GC statistics
        """
        gc_profile = self.calculate_gc_profile(sequence)
        overall_gc = self.calculate_gc_content(sequence)
        gc_variance = self.calculate_gc_variance(sequence)
        violations = self.find_gc_violations(sequence)
        gc_score = self.calculate_gc_score(sequence)
        
        return {
            'overall_gc': overall_gc,
            'mean_gc': float(np.mean(gc_profile)) if gc_profile else 0.0,
            'min_gc': float(np.min(gc_profile)) if gc_profile else 0.0,
            'max_gc': float(np.max(gc_profile)) if gc_profile else 0.0,
            'std_gc': float(np.std(gc_profile)) if gc_profile else 0.0,
            'variance_gc': gc_variance,
            'violation_count': len(violations),
            'violation_fraction': len(violations) / len(gc_profile) if gc_profile else 0.0,
            'gc_score': gc_score
        }
    
    def check_extreme_regions(self, sequence: str) -> List[Tuple[int, int, str]]:
        """
        Check for extended regions with extreme GC content.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            List of (start, end, type) where type is "low_gc" or "high_gc"
        """
        gc_profile = self.calculate_gc_profile(sequence)
        extreme_regions = []
        
        current_region_start = None
        current_region_type = None
        
        for i, gc in enumerate(gc_profile):
            if gc < self.gc_min:
                if current_region_type != "low_gc":
                    if current_region_start is not None:
                        # End previous region
                        extreme_regions.append((
                            current_region_start,
                            i + self.window_size - 1,
                            current_region_type
                        ))
                    current_region_start = i
                    current_region_type = "low_gc"
            elif gc > self.gc_max:
                if current_region_type != "high_gc":
                    if current_region_start is not None:
                        # End previous region
                        extreme_regions.append((
                            current_region_start,
                            i + self.window_size - 1,
                            current_region_type
                        ))
                    current_region_start = i
                    current_region_type = "high_gc"
            else:
                if current_region_start is not None:
                    # End current region
                    extreme_regions.append((
                        current_region_start,
                        i + self.window_size - 1,
                        current_region_type
                    ))
                    current_region_start = None
                    current_region_type = None
        
        # Handle region extending to end
        if current_region_start is not None:
            extreme_regions.append((
                current_region_start,
                len(sequence) - 1,
                current_region_type
            ))
        
        return extreme_regions




