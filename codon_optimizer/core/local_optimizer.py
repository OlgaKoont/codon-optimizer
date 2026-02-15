"""
Local optimization approach inspired by DnaChisel.
Instead of generating all combinations, we iteratively optimize problematic regions.
"""

from typing import List, Tuple, Dict, Optional, Set
import random
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.core.scorer import MultiCriteriaScorer
from codon_optimizer.core.codon_generator import CodonCombinationGenerator
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from codon_optimizer.analysis.gc_content import GCContentAnalyzer
from codon_optimizer.analysis.motifs import MotifDetector


class LocalOptimizationRegion:
    """Represents a region that needs local optimization."""
    
    def __init__(self, start: int, end: int, reason: str, priority: float = 1.0):
        """
        Initialize optimization region.
        
        Args:
            start: Start position (codon index)
            end: End position (codon index, exclusive)
            reason: Reason for optimization (e.g., "low_cai", "gc_violation", "motif")
            priority: Priority of this region (higher = more important)
        """
        self.start = start
        self.end = end
        self.reason = reason
        self.priority = priority
    
    def __repr__(self):
        return f"Region({self.start}-{self.end}, {self.reason}, priority={self.priority})"


class LocalCodonOptimizer:
    """
    Local optimizer inspired by DnaChisel's approach.
    
    Instead of generating all combinations, this optimizer:
    1. Starts with an optimal sequence
    2. Identifies problematic regions
    3. Locally optimizes each region
    4. Iteratively improves until convergence
    """
    
    def __init__(self, config: OptimizationConfig):
        """
        Initialize local optimizer.
        
        Args:
            config: Optimization configuration
        """
        self.config = config
        self.scorer = MultiCriteriaScorer(config)
        self.codon_generator = CodonCombinationGenerator(host=config.host)
        self.codon_analyzer = CodonUsageAnalyzer(
            host=config.host,
            max_codon_usage=config.max_codon_usage
        )
        self.gc_analyzer = GCContentAnalyzer(
            gc_min=config.gc_min,
            gc_max=config.gc_max,
            window_size=config.gc_window_size
        )
        self.motif_detector = MotifDetector()
    
    def _generate_initial_sequence(self, protein_sequence: str) -> str:
        """
        Generate initial sequence using optimal codons.
        
        Args:
            protein_sequence: Protein sequence
            
        Returns:
            Initial DNA sequence
        """
        dna_sequence = []
        for aa in protein_sequence:
            optimal_codon = self.codon_analyzer.get_optimal_codon(aa)
            dna_sequence.append(optimal_codon)
        return ''.join(dna_sequence)
    
    def _identify_problematic_regions(self, 
                                     dna_sequence: str,
                                     protein_sequence: str,
                                     window_size: int = 10) -> List[LocalOptimizationRegion]:
        """
        Identify regions that need optimization.
        
        Args:
            dna_sequence: Current DNA sequence
            protein_sequence: Protein sequence
            window_size: Size of window for analysis (in codons)
            
        Returns:
            List of regions that need optimization, sorted by priority
        """
        regions = []
        num_codons = len(protein_sequence)
        
        # Analyze in windows
        for i in range(0, num_codons, window_size // 2):  # Overlapping windows
            end = min(i + window_size, num_codons)
            if end <= i:
                continue
            
            # Extract region
            region_dna = dna_sequence[i*3:(end*3)]
            region_protein = protein_sequence[i:end]
            
            if len(region_dna) < 3:
                continue
            
            # Check CAI
            cai = self.codon_analyzer.calculate_cai(region_dna)
            if cai < self.config.target_cai_min:
                priority = (self.config.target_cai_min - cai) * 2.0
                regions.append(LocalOptimizationRegion(
                    i, end, f"low_cai_{cai:.3f}", priority
                ))
            
            # Check GC content
            gc_stats = self.gc_analyzer.analyze_gc_content(region_dna)
            overall_gc = gc_stats['overall_gc']
            if overall_gc < self.config.gc_min or overall_gc > self.config.gc_max:
                priority = abs(overall_gc - (self.config.gc_min + self.config.gc_max) / 2) * 1.5
                regions.append(LocalOptimizationRegion(
                    i, end, f"gc_violation_{overall_gc:.3f}", priority
                ))
            
            # Check for problematic motifs
            all_motifs = self.motif_detector.find_motifs(region_dna)
            total_motifs = sum(len(motif_list) for motif_list in all_motifs.values())
            if total_motifs > 0:
                priority = total_motifs * 0.5
                regions.append(LocalOptimizationRegion(
                    i, end, f"motifs_{total_motifs}", priority
                ))
        
        # Sort by priority (highest first)
        regions.sort(key=lambda r: r.priority, reverse=True)
        
        # Merge overlapping regions
        merged_regions = self._merge_overlapping_regions(regions)
        
        return merged_regions
    
    def _merge_overlapping_regions(self, regions: List[LocalOptimizationRegion]) -> List[LocalOptimizationRegion]:
        """
        Merge overlapping regions.
        
        Args:
            regions: List of regions
            
        Returns:
            Merged regions
        """
        if not regions:
            return []
        
        merged = []
        current = regions[0]
        
        for region in regions[1:]:
            if region.start <= current.end:
                # Overlapping - merge
                current.end = max(current.end, region.end)
                current.priority = max(current.priority, region.priority)
                current.reason = f"{current.reason}+{region.reason}"
            else:
                merged.append(current)
                current = region
        
        merged.append(current)
        return merged
    
    def _optimize_region(self,
                        dna_sequence: str,
                        protein_sequence: str,
                        region: LocalOptimizationRegion,
                        max_combinations: int = 100) -> Tuple[str, float]:
        """
        Optimize a specific region locally.
        
        Args:
            dna_sequence: Current full DNA sequence
            protein_sequence: Full protein sequence
            region: Region to optimize
            max_combinations: Maximum combinations to try for this region
            
        Returns:
            Tuple of (optimized_sequence, improvement_score)
        """
        # Extract region sequences
        region_protein = protein_sequence[region.start:region.end]
        region_dna_start = region.start * 3
        region_dna_end = region.end * 3
        region_dna = dna_sequence[region_dna_start:region_dna_end]
        
        # Generate codon options for this region
        codon_options_list = []
        for aa in region_protein:
            codons = self.codon_generator.get_all_codons_for_aa(aa)
            codon_options_list.append(codons)
        
        # Limit combinations if too many
        total_region_combinations = 1
        for codons in codon_options_list:
            total_region_combinations *= len(codons)
        
        if total_region_combinations > max_combinations:
            # Use smart generation
            combinations = list(self.codon_generator.generate_smart_combinations(
                region_protein,
                max_combinations=max_combinations,
                strategy='diverse'
            ))
        else:
            # Generate all
            combinations = list(self.codon_generator.generate_all_combinations(region_protein))
        
        # Evaluate each combination in context
        best_region_dna = region_dna
        best_score = -float('inf')
        
        for candidate_region_dna in combinations:
            # Reconstruct full sequence
            candidate_full = (
                dna_sequence[:region_dna_start] +
                candidate_region_dna +
                dna_sequence[region_dna_end:]
            )
            
            # Score this candidate
            score_dict = self.scorer.calculate_total_score(candidate_full, protein_sequence)
            score = score_dict['total_score']
            
            if score > best_score:
                best_score = score
                best_region_dna = candidate_region_dna
        
        # Calculate improvement
        original_score_dict = self.scorer.calculate_total_score(dna_sequence, protein_sequence)
        original_score = original_score_dict['total_score']
        improvement = best_score - original_score
        
        # Reconstruct optimized sequence
        optimized_sequence = (
            dna_sequence[:region_dna_start] +
            best_region_dna +
            dna_sequence[region_dna_end:]
        )
        
        return optimized_sequence, improvement
    
    def optimize(self,
                protein_sequence: str,
                max_iterations: int = 50,
                convergence_threshold: float = 0.001,
                max_region_combinations: int = 100) -> Tuple[str, Dict]:
        """
        Optimize sequence using local iterative approach.
        
        Args:
            protein_sequence: Protein sequence to optimize
            max_iterations: Maximum number of optimization iterations
            convergence_threshold: Minimum improvement to continue
            max_region_combinations: Max combinations to try per region
            
        Returns:
            Tuple of (optimized_dna_sequence, metadata_dict)
        """
        # Start with optimal sequence
        current_sequence = self._generate_initial_sequence(protein_sequence)
        current_score = self.scorer.calculate_total_score(current_sequence, protein_sequence)['total_score']
        
        iteration = 0
        total_improvements = 0
        regions_optimized = 0
        
        while iteration < max_iterations:
            iteration += 1
            
            # Identify problematic regions
            regions = self._identify_problematic_regions(current_sequence, protein_sequence)
            
            if not regions:
                # No problematic regions - we're done!
                break
            
            # Optimize top regions (up to 5 per iteration)
            improved = False
            for region in regions[:5]:  # Limit to top 5 regions per iteration
                optimized_sequence, improvement = self._optimize_region(
                    current_sequence,
                    protein_sequence,
                    region,
                    max_combinations=max_region_combinations
                )
                
                if improvement > convergence_threshold:
                    current_sequence = optimized_sequence
                    current_score += improvement
                    total_improvements += improvement
                    regions_optimized += 1
                    improved = True
            
            # Check convergence
            if not improved:
                # No improvement in this iteration
                break
        
        # Final score
        final_score_dict = self.scorer.calculate_total_score(current_sequence, protein_sequence)
        final_score = final_score_dict['total_score']
        
        metadata = {
            'iterations': iteration,
            'regions_optimized': regions_optimized,
            'total_improvements': total_improvements,
            'final_score': final_score,
            'initial_score': current_score - total_improvements,
            'optimization_method': 'local_iterative'
        }
        
        return current_sequence, metadata

