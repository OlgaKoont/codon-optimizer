"""
Protein integrity constraints for codon optimization.
"""

from typing import List, Tuple, Dict, Set
from Bio.Seq import Seq
from Bio.SeqUtils import seq1


class ProteinIntegrityChecker:
    """Checker for protein integrity constraints."""
    
    def __init__(self,
                 protect_cysteines: bool = True,
                 protect_glycosylation: bool = True,
                 protect_cdr_regions: bool = True):
        """
        Initialize protein integrity checker.
        
        Args:
            protect_cysteines: Protect canonical cysteines from changes
            protect_glycosylation: Protect N-glycosylation sites
            protect_cdr_regions: Protect CDR regions in antibodies
        """
        self.protect_cysteines = protect_cysteines
        self.protect_glycosylation = protect_glycosylation
        self.protect_cdr_regions = protect_cdr_regions
        
        # N-glycosylation motif: NXS/T where X is not P
        self.nglycosylation_pattern = r'N[^P][ST]'
    
    def validate_protein_integrity(self,
                                  original_protein: str,
                                  optimized_protein: str) -> Tuple[bool, List[str]]:
        """
        Validate that optimized protein matches original.
        
        Args:
            original_protein: Original protein sequence
            optimized_protein: Optimized protein sequence
            
        Returns:
            Tuple of (is_valid, list_of_errors)
        """
        errors = []
        
        if len(original_protein) != len(optimized_protein):
            errors.append(f"Length mismatch: {len(original_protein)} vs {len(optimized_protein)}")
            return False, errors
        
        mismatches = []
        for i, (orig_aa, opt_aa) in enumerate(zip(original_protein, optimized_protein)):
            if orig_aa != opt_aa:
                mismatches.append(f"Position {i}: {orig_aa} -> {opt_aa}")
        
        if mismatches:
            errors.extend(mismatches)
            return False, errors
        
        return True, []
    
    def find_cysteines(self, protein_sequence: str) -> List[int]:
        """
        Find positions of cysteines in protein.
        
        Args:
            protein_sequence: Protein sequence
            
        Returns:
            List of cysteine positions
        """
        return [i for i, aa in enumerate(protein_sequence) if aa == 'C']
    
    def find_glycosylation_sites(self, protein_sequence: str) -> List[Tuple[int, str]]:
        """
        Find N-glycosylation sites (NXS/T) in protein.
        
        Args:
            protein_sequence: Protein sequence
            
        Returns:
            List of (position, motif) tuples
        """
        import re
        sites = []
        
        for match in re.finditer(self.nglycosylation_pattern, protein_sequence):
            position = match.start()
            motif = match.group()
            sites.append((position, motif))
        
        return sites
    
    def find_unpaired_cysteines(self, protein_sequence: str) -> List[int]:
        """
        Find potentially unpaired cysteines.
        This is a simplified check - real analysis would require structure.
        
        Args:
            protein_sequence: Protein sequence
            
        Returns:
            List of positions of potentially unpaired cysteines
        """
        cysteines = self.find_cysteines(protein_sequence)
        
        # Simple heuristic: odd number of cysteines suggests unpaired
        # In practice, this would use structure prediction
        if len(cysteines) % 2 != 0:
            return cysteines
        else:
            return []
    
    def find_hydrophobic_regions(self,
                                 protein_sequence: str,
                                 min_length: int = 10) -> List[Tuple[int, int]]:
        """
        Find long hydrophobic regions.
        
        Args:
            protein_sequence: Protein sequence
            min_length: Minimum length of hydrophobic region
            
        Returns:
            List of (start, end) tuples
        """
        hydrophobic_aas = set('AILMVFYW')
        regions = []
        
        i = 0
        while i < len(protein_sequence):
            if protein_sequence[i] in hydrophobic_aas:
                start = i
                while i < len(protein_sequence) and protein_sequence[i] in hydrophobic_aas:
                    i += 1
                length = i - start
                if length >= min_length:
                    regions.append((start, i))
            else:
                i += 1
        
        return regions
    
    def get_protected_positions(self,
                               protein_sequence: str,
                               cdr_regions: List[Tuple[int, int]] = None) -> Set[int]:
        """
        Get set of positions that must be protected from codon changes.
        
        Args:
            protein_sequence: Protein sequence
            cdr_regions: List of (start, end) tuples for CDR regions (for antibodies)
            
        Returns:
            Set of protected positions
        """
        protected = set()
        
        # Protect cysteines
        if self.protect_cysteines:
            protected.update(self.find_cysteines(protein_sequence))
        
        # Protect glycosylation sites
        if self.protect_glycosylation:
            for pos, _ in self.find_glycosylation_sites(protein_sequence):
                # Protect all 3 positions of motif
                protected.update(range(pos, pos + 3))
        
        # Protect CDR regions
        if self.protect_cdr_regions and cdr_regions:
            for start, end in cdr_regions:
                protected.update(range(start, end))
        
        return protected
    
    def check_integrity_violations(self,
                                  original_protein: str,
                                  optimized_protein: str,
                                  cdr_regions: List[Tuple[int, int]] = None) -> Dict[str, any]:
        """
        Check for protein integrity violations.
        
        Args:
            original_protein: Original protein sequence
            optimized_protein: Optimized protein sequence
            cdr_regions: CDR regions for antibodies
            
        Returns:
            Dictionary with violation information
        """
        is_valid, errors = self.validate_protein_integrity(original_protein, optimized_protein)
        
        violations = {
            'is_valid': is_valid,
            'errors': errors,
            'cysteine_changes': [],
            'glycosylation_changes': [],
            'cdr_changes': [],
        }
        
        if not is_valid:
            # Check specific violations
            orig_cys = set(self.find_cysteines(original_protein))
            opt_cys = set(self.find_cysteines(optimized_protein))
            
            if orig_cys != opt_cys:
                violations['cysteine_changes'] = [
                    f"Lost: {orig_cys - opt_cys}, Gained: {opt_cys - orig_cys}"
                ]
            
            orig_glyc = set(self.find_glycosylation_sites(original_protein))
            opt_glyc = set(self.find_glycosylation_sites(optimized_protein))
            
            if orig_glyc != opt_glyc:
                violations['glycosylation_changes'] = [
                    f"Lost: {orig_glyc - opt_glyc}, Gained: {opt_glyc - orig_glyc}"
                ]
        
        return violations
    
    def mark_problematic_motifs(self, protein_sequence: str) -> Dict[str, List]:
        """
        Mark motifs associated with poor physicochemical properties.
        
        Args:
            protein_sequence: Protein sequence
            
        Returns:
            Dictionary with marked motifs
        """
        return {
            'unpaired_cysteines': self.find_unpaired_cysteines(protein_sequence),
            'hydrophobic_regions': self.find_hydrophobic_regions(protein_sequence),
            'glycosylation_sites': self.find_glycosylation_sites(protein_sequence),
        }




