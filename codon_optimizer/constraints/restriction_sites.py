"""
Restriction enzyme site detection and management.
"""

from typing import List, Tuple, Dict, Set
from Bio.Restriction import Restriction, RestrictionBatch
from Bio.Seq import Seq


# Common restriction enzymes used in cloning
COMMON_RESTRICTION_ENZYMES = {
    'BsaI': 'GGTCTC',
    'BsmBI': 'CGTCTC',
    'NotI': 'GCGGCCGC',
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC',
    'HindIII': 'AAGCTT',
    'XhoI': 'CTCGAG',
    'SalI': 'GTCGAC',
    'NcoI': 'CCATGG',
    'AgeI': 'ACCGGT',
    'KpnI': 'GGTACC',
    'XbaI': 'TCTAGA',
    'SpeI': 'ACTAGT',
    'PstI': 'CTGCAG',
    'SphI': 'GCATGC',
    'SacI': 'GAGCTC',
    'SacII': 'CCGCGG',
    'SmaI': 'CCCGGG',
    'XmaI': 'CCCGGG',
    'ApaI': 'GGGCCC',
    'SwaI': 'ATTTAAAT',
    'PacI': 'TTAATTAA',
}


class RestrictionSiteManager:
    """Manager for restriction enzyme sites."""
    
    def __init__(self,
                 sites_to_remove: List[str] = None,
                 sites_to_keep: List[str] = None):
        """
        Initialize restriction site manager.
        
        Args:
            sites_to_remove: List of enzyme names to remove from sequence
            sites_to_keep: List of enzyme names that must be preserved
        """
        self.sites_to_remove = set(sites_to_remove or [])
        self.sites_to_keep = set(sites_to_keep or [])
        
        # Build recognition sequences
        self.recognition_sequences = {}
        for enzyme_name, sequence in COMMON_RESTRICTION_ENZYMES.items():
            self.recognition_sequences[enzyme_name] = sequence
    
    def find_restriction_sites(self, sequence: str, enzyme_name: str = None) -> Dict[str, List[int]]:
        """
        Find restriction sites in sequence.
        
        Args:
            sequence: DNA sequence
            enzyme_name: Specific enzyme to search (None for all)
            
        Returns:
            Dictionary mapping enzyme names to list of positions
        """
        sequence = sequence.upper()
        results = {}
        
        enzymes_to_search = [enzyme_name] if enzyme_name else self.recognition_sequences.keys()
        
        for enzyme in enzymes_to_search:
            if enzyme not in self.recognition_sequences:
                continue
            
            recognition_seq = self.recognition_sequences[enzyme]
            positions = []
            
            # Find all occurrences
            start = 0
            while True:
                pos = sequence.find(recognition_seq, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1
            
            if positions:
                results[enzyme] = positions
        
        return results
    
    def find_sites_to_remove(self, sequence: str) -> Dict[str, List[int]]:
        """Find restriction sites that should be removed."""
        all_sites = self.find_restriction_sites(sequence)
        return {enzyme: positions 
                for enzyme, positions in all_sites.items() 
                if enzyme in self.sites_to_remove}
    
    def find_sites_to_keep(self, sequence: str) -> Dict[str, List[int]]:
        """Find restriction sites that must be preserved."""
        all_sites = self.find_restriction_sites(sequence)
        return {enzyme: positions 
                for enzyme, positions in all_sites.items() 
                if enzyme in self.sites_to_keep}
    
    def check_violations(self, sequence: str) -> Dict[str, any]:
        """
        Check for restriction site violations.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with violation information
        """
        sites_to_remove = self.find_sites_to_remove(sequence)
        sites_to_keep = self.find_sites_to_keep(sequence)
        
        # Check if sites to keep are present
        missing_required = []
        for enzyme in self.sites_to_keep:
            if enzyme not in sites_to_keep or not sites_to_keep[enzyme]:
                missing_required.append(enzyme)
        
        return {
            'sites_to_remove': sites_to_remove,
            'sites_to_keep': sites_to_keep,
            'has_unwanted_sites': len(sites_to_remove) > 0,
            'missing_required_sites': missing_required,
            'violation_count': sum(len(positions) for positions in sites_to_remove.values())
        }
    
    def calculate_cloning_score(self, sequence: str) -> float:
        """
        Calculate cloning compatibility score (0-1, higher is better).
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Score between 0 and 1
        """
        violations = self.check_violations(sequence)
        
        score = 1.0
        
        # Penalty for unwanted sites
        unwanted_count = violations['violation_count']
        if unwanted_count > 0:
            # Penalty increases with number of unwanted sites
            penalty = min(0.5, unwanted_count * 0.1)
            score -= penalty
        
        # Penalty for missing required sites
        missing_count = len(violations['missing_required_sites'])
        if missing_count > 0:
            penalty = min(0.5, missing_count * 0.2)
            score -= penalty
        
        return max(0.0, score)
    
    def suggest_codon_changes(self,
                             sequence: str,
                             position: int,
                             amino_acid: str,
                             codon_usage_analyzer) -> List[str]:
        """
        Suggest codon changes to remove restriction site.
        
        Args:
            sequence: DNA sequence
            position: Position of codon to change
            amino_acid: Amino acid at this position
            codon_usage_analyzer: CodonUsageAnalyzer instance
            
        Returns:
            List of alternative codons that avoid restriction sites
        """
        # Get alternative codons
        alternatives = codon_usage_analyzer.get_codon_options(amino_acid)
        
        # Check which alternatives avoid restriction sites
        safe_codons = []
        
        for alt_codon in alternatives:
            # Create test sequence with this codon
            test_seq = sequence[:position] + alt_codon + sequence[position + 3:]
            
            # Check if this creates unwanted sites
            unwanted = self.find_sites_to_remove(test_seq)
            if not unwanted:
                safe_codons.append(alt_codon)
        
        return safe_codons if safe_codons else alternatives




