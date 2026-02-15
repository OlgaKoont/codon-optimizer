"""
mRNA secondary structure analysis using ViennaRNA.
"""

from typing import Tuple, Optional, Dict, List
import subprocess
import tempfile
import os
from pathlib import Path


class MRNAStructureAnalyzer:
    """Analyzer for mRNA secondary structure using ViennaRNA."""
    
    def __init__(self,
                 max_5prime_dg: float = -10.0,
                 structure_region_start: int = -50,
                 structure_region_end: int = 100):
        """
        Initialize mRNA structure analyzer.
        
        Args:
            max_5prime_dg: Maximum allowed ΔG for 5' region (kcal/mol)
            structure_region_start: Start of critical region relative to start codon
            structure_region_end: End of critical region relative to start codon
        """
        self.max_5prime_dg = max_5prime_dg
        self.structure_region_start = structure_region_start
        self.structure_region_end = structure_region_end
        self.rnafold_available = self._check_rnafold_available()
    
    def _check_rnafold_available(self) -> bool:
        """Check if RNAfold is available in PATH."""
        try:
            result = subprocess.run(
                ['RNAfold', '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
    
    def _call_rnafold(self, sequence: str, temperature: float = 37.0) -> Tuple[Optional[float], Optional[str]]:
        """
        Call RNAfold to calculate secondary structure and ΔG.
        
        Args:
            sequence: RNA sequence
            temperature: Temperature in Celsius
            
        Returns:
            Tuple of (ΔG in kcal/mol, dot-bracket structure) or (None, None) if failed
        """
        if not self.rnafold_available:
            # Fallback: estimate ΔG based on sequence length
            # This is a rough approximation
            estimated_dg = -0.1 * len(sequence)  # Rough estimate
            return estimated_dg, None
        
        try:
            # Convert DNA to RNA if needed
            rna_seq = sequence.replace('T', 'U').upper()
            
            # Create temporary file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
                f.write(f">temp\n{rna_seq}\n")
                temp_file = f.name
            
            try:
                # Call RNAfold
                result = subprocess.run(
                    ['RNAfold', '--temp', str(temperature), '--noPS', temp_file],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
                
                if result.returncode == 0:
                    # Parse output
                    lines = result.stdout.strip().split('\n')
                    if len(lines) >= 2:
                        # Second line contains sequence and structure
                        structure_line = lines[1]
                        # Extract ΔG from last part (format: "sequence (structure) dG_value")
                        parts = structure_line.split()
                        if len(parts) >= 2:
                            dg_str = parts[-1].strip('()')
                            try:
                                dg = float(dg_str)
                                # Extract structure (between parentheses)
                                structure = None
                                if '(' in structure_line and ')' in structure_line:
                                    start = structure_line.find('(')
                                    end = structure_line.rfind(')')
                                    if start < end:
                                        structure = structure_line[start+1:end]
                                
                                return dg, structure
                            except ValueError:
                                pass
                
                return None, None
                
            finally:
                # Clean up temp file
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
                    
        except Exception as e:
            print(f"Error calling RNAfold: {e}")
            return None, None
    
    def calculate_structure_dg(self, 
                              sequence: str,
                              region_start: Optional[int] = None,
                              region_end: Optional[int] = None) -> Optional[float]:
        """
        Calculate ΔG for a specific region of the sequence.
        
        Args:
            sequence: DNA sequence
            region_start: Start position (None for full sequence)
            region_end: End position (None for full sequence)
            
        Returns:
            ΔG in kcal/mol or None if calculation failed
        """
        if region_start is not None and region_end is not None:
            region_seq = sequence[region_start:region_end]
        else:
            region_seq = sequence
        
        if not region_seq:
            return None
        
        dg, _ = self._call_rnafold(region_seq)
        return dg
    
    def analyze_5prime_region(self, 
                             sequence: str,
                             start_codon_pos: int = 0) -> Dict[str, Optional[float]]:
        """
        Analyze 5' region for stable secondary structures.
        
        Args:
            sequence: DNA sequence
            start_codon_pos: Position of start codon (0-indexed)
            
        Returns:
            Dictionary with structure analysis results
        """
        # Calculate region boundaries
        region_start = max(0, start_codon_pos + self.structure_region_start)
        region_end = min(len(sequence), start_codon_pos + self.structure_region_end)
        
        # Calculate ΔG for critical region
        critical_dg = self.calculate_structure_dg(sequence, region_start, region_end)
        
        # Calculate ΔG for full sequence
        full_dg = self.calculate_structure_dg(sequence)
        
        # Calculate ΔG for 5' UTR-like region (before start codon)
        if start_codon_pos > 0:
            utr_start = max(0, start_codon_pos - 50)
            utr_dg = self.calculate_structure_dg(sequence, utr_start, start_codon_pos)
        else:
            utr_dg = None
        
        return {
            'critical_region_dg': critical_dg,
            'full_sequence_dg': full_dg,
            'utr_dg': utr_dg,
            'critical_region_start': region_start,
            'critical_region_end': region_end,
            'has_stable_structure': critical_dg is not None and critical_dg < self.max_5prime_dg
        }
    
    def calculate_structure_score(self, sequence: str, start_codon_pos: int = 0) -> float:
        """
        Calculate structure score (0-1, higher is better).
        
        Args:
            sequence: DNA sequence
            start_codon_pos: Position of start codon
            
        Returns:
            Score between 0 and 1
        """
        analysis = self.analyze_5prime_region(sequence, start_codon_pos)
        critical_dg = analysis['critical_region_dg']
        
        if critical_dg is None:
            return 0.5  # Unknown structure
        
        # Score based on how far ΔG is from threshold
        if critical_dg >= self.max_5prime_dg:
            # Good: structure is not too stable
            score = 1.0
        else:
            # Penalty for stable structures
            # More negative ΔG = more stable = worse
            penalty = (self.max_5prime_dg - critical_dg) / abs(self.max_5prime_dg)
            score = max(0.0, 1.0 - penalty)
        
        return score
    
    def find_stable_structures(self, 
                              sequence: str,
                              window_size: int = 50,
                              step_size: int = 10) -> List[Tuple[int, int, float]]:
        """
        Find regions with stable secondary structures.
        
        Args:
            sequence: DNA sequence
            window_size: Size of sliding window
            step_size: Step size for sliding window
            
        Returns:
            List of (start, end, dg) tuples for stable structures
        """
        stable_regions = []
        
        for i in range(0, len(sequence) - window_size, step_size):
            window_seq = sequence[i:i + window_size]
            dg = self.calculate_structure_dg(window_seq)
            
            if dg is not None and dg < self.max_5prime_dg:
                stable_regions.append((i, i + window_size, dg))
        
        return stable_regions
    
    def get_structure_statistics(self, sequence: str, start_codon_pos: int = 0) -> Dict:
        """
        Get comprehensive structure statistics.
        
        Args:
            sequence: DNA sequence
            start_codon_pos: Position of start codon
            
        Returns:
            Dictionary with structure statistics
        """
        analysis = self.analyze_5prime_region(sequence, start_codon_pos)
        structure_score = self.calculate_structure_score(sequence, start_codon_pos)
        stable_structures = self.find_stable_structures(sequence)
        
        return {
            'critical_region_dg': analysis['critical_region_dg'],
            'full_sequence_dg': analysis['full_sequence_dg'],
            'utr_dg': analysis['utr_dg'],
            'has_stable_structure': analysis['has_stable_structure'],
            'structure_score': structure_score,
            'stable_structure_count': len(stable_structures),
            'stable_structures': stable_structures
        }
    
    def suggest_sequence_changes(self,
                                sequence: str,
                                position: int,
                                target_dg: float) -> List[str]:
        """
        Suggest sequence changes to reduce structure stability.
        This is a placeholder - actual implementation would require
        integration with the optimizer.
        
        Args:
            sequence: DNA sequence
            position: Position to modify
            target_dg: Target ΔG value
            
        Returns:
            List of suggested alternative sequences
        """
        # This will be used by optimizer
        # For now, return empty list
        return []




