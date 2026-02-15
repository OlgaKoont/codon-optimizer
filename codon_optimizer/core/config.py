"""
Configuration settings for codon optimization.
"""

from dataclasses import dataclass
from typing import List, Optional, Dict


@dataclass
class OptimizationConfig:
    """Configuration for codon optimization."""
    
    # Host organism
    host: str = "CHO"
    
    # CAI settings
    target_cai_min: float = 0.7
    target_cai_max: float = 0.9
    max_codon_usage: float = 0.5  # Maximum usage of single codon to avoid tRNA depletion
    
    # GC content settings
    gc_min: float = 0.30
    gc_max: float = 0.80
    gc_window_size: int = 50
    
    # mRNA structure settings
    max_5prime_dg: float = -10.0  # Maximum ΔG for 5' region (kcal/mol)
    structure_region_start: int = -50  # Start of critical region relative to start codon
    structure_region_end: int = 100  # End of critical region relative to start codon
    
    # Expression settings
    target_expression_min: float = 0.7
    target_expression_max: float = 0.9
    
    # Restriction sites
    restriction_sites_to_remove: List[str] = None
    restriction_sites_to_keep: List[str] = None
    
    # Optimization algorithm
    algorithm: str = "genetic"  # "genetic" or "simulated_annealing"
    population_size: int = 50
    max_generations: int = 100
    mutation_rate: float = 0.1
    
    # Pareto optimization
    use_pareto: bool = True
    pareto_solutions: int = 10
    
    # Antibody-specific
    is_antibody: bool = False
    heavy_chain_id: Optional[str] = None
    light_chain_id: Optional[str] = None
    balance_chains: bool = True
    
    # Scoring weights
    weight_cai: float = 0.25
    weight_gc: float = 0.15
    weight_codon_pairs: float = 0.15
    weight_mrna_structure: float = 0.20
    weight_motifs: float = 0.15
    weight_cloning: float = 0.10
    
    # Kozak sequence
    kozak_sequence: str = "GCCACC"
    
    def __post_init__(self):
        """Initialize default values for lists."""
        if self.restriction_sites_to_remove is None:
            self.restriction_sites_to_remove = []
        if self.restriction_sites_to_keep is None:
            self.restriction_sites_to_keep = []
    
    def validate(self):
        """Validate configuration parameters."""
        if not 0 <= self.target_cai_min <= self.target_cai_max <= 1.0:
            raise ValueError("CAI values must be between 0 and 1")
        
        if not 0 <= self.gc_min <= self.gc_max <= 1.0:
            raise ValueError("GC content values must be between 0 and 1")
        
        if sum([self.weight_cai, self.weight_gc, self.weight_codon_pairs,
                self.weight_mrna_structure, self.weight_motifs, self.weight_cloning]) != 1.0:
            raise ValueError("Weights must sum to 1.0")
        
        if self.population_size < 10:
            raise ValueError("Population size must be at least 10")
        
        if self.max_generations < 1:
            raise ValueError("Max generations must be at least 1")




