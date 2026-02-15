"""
Detection of problematic sequence motifs.
"""

import re
from typing import List, Tuple, Dict, Set
from collections import defaultdict


class MotifDetector:
    """Detector for problematic sequence motifs."""
    
    def __init__(self):
        """Initialize motif detector with pattern definitions."""
        self.motif_patterns = self._build_motif_patterns()
    
    def _build_motif_patterns(self) -> Dict[str, List[str]]:
        """Build dictionary of motif patterns."""
        return {
            # TATA boxes
            'tata_box': [
                r'TATA[AT][AT]A',  # Standard TATA box
                r'TATAAA',
                r'TATATA',
            ],
            
            # Cryptic promoters (simplified - real promoters are more complex)
            'cryptic_promoter': [
                r'[AT]TATA[AT][AT]',  # TATA variants
                r'CAAT',  # CAAT box
                r'GGGCGG',  # GC box
            ],
            
            # Chi sites (E.coli, adapted for CHO)
            'chi_site': [
                r'GCTGGTGG',  # E.coli chi site
                r'[GC]CTGG[TC]GG',  # Variants
            ],
            
            # Internal RBS (Shine-Dalgarno-like)
            'internal_rbs': [
                r'AGGAGG',  # Strong SD
                r'GGAGG',   # Medium SD
                r'AGGA',    # Weak SD
            ],
            
            # ARE (AU-rich elements) - for RNA
            'are_motif': [
                r'ATTTA',   # Classic ARE
                r'ATTT[AT]A',
                r'[AT]TTTA',
            ],
            
            # Splice sites
            'splice_donor': [
                r'GT[AG]',  # Donor site
            ],
            'splice_acceptor': [
                r'[CT]AG',  # Acceptor site
            ],
            'branch_point': [
                r'[CT]T[AG]AC',  # Branch point consensus
                r'YTRAY',  # Y = pyrimidine, R = purine
            ],
            
            # Polyadenylation signals
            'polyA_signal': [
                r'AATAAA',  # Canonical
                r'ATTAAA',  # Variant
                r'AATACA',  # Variant
                r'ACTAAA',  # Variant
            ],
        }
    
    def find_motifs(self, sequence: str, motif_type: str = None) -> Dict[str, List[Tuple[int, int, str]]]:
        """
        Find all motifs in sequence.
        
        Args:
            sequence: DNA sequence
            motif_type: Specific motif type to search (None for all)
            
        Returns:
            Dictionary mapping motif type to list of (start, end, match) tuples
        """
        sequence = sequence.upper()
        results = defaultdict(list)
        
        patterns_to_search = {}
        if motif_type:
            if motif_type in self.motif_patterns:
                patterns_to_search[motif_type] = self.motif_patterns[motif_type]
        else:
            patterns_to_search = self.motif_patterns
        
        for motif_type, patterns in patterns_to_search.items():
            for pattern in patterns:
                for match in re.finditer(pattern, sequence, re.IGNORECASE):
                    start = match.start()
                    end = match.end()
                    match_str = match.group()
                    results[motif_type].append((start, end, match_str))
        
        return dict(results)
    
    def find_tata_boxes(self, sequence: str) -> List[Tuple[int, int, str]]:
        """Find TATA boxes in sequence."""
        return self.find_motifs(sequence, 'tata_box').get('tata_box', [])
    
    def find_cryptic_promoters(self, sequence: str) -> List[Tuple[int, int, str]]:
        """Find cryptic promoters in sequence."""
        return self.find_motifs(sequence, 'cryptic_promoter').get('cryptic_promoter', [])
    
    def find_splice_sites(self, sequence: str) -> Dict[str, List[Tuple[int, int, str]]]:
        """Find splice sites in sequence."""
        motifs = self.find_motifs(sequence)
        return {
            'donor': motifs.get('splice_donor', []),
            'acceptor': motifs.get('splice_acceptor', []),
            'branch_point': motifs.get('branch_point', [])
        }
    
    def find_polyA_signals(self, sequence: str) -> List[Tuple[int, int, str]]:
        """Find polyadenylation signals in sequence."""
        return self.find_motifs(sequence, 'polyA_signal').get('polyA_signal', [])
    
    def find_homopolymers(self, sequence: str, min_length: int = 5) -> List[Tuple[int, int, str]]:
        """
        Find homopolymer tracts (repeats of same nucleotide).
        
        Args:
            sequence: DNA sequence
            min_length: Minimum length of homopolymer
            
        Returns:
            List of (start, end, nucleotide) tuples
        """
        sequence = sequence.upper()
        homopolymers = []
        
        i = 0
        while i < len(sequence):
            nucleotide = sequence[i]
            start = i
            length = 1
            
            # Count consecutive same nucleotides
            while i + length < len(sequence) and sequence[i + length] == nucleotide:
                length += 1
            
            if length >= min_length:
                homopolymers.append((start, start + length, nucleotide * length))
            
            i += length
        
        return homopolymers
    
    def find_tandem_repeats(self, 
                           sequence: str,
                           min_repeat_length: int = 3,
                           min_repeats: int = 3) -> List[Tuple[int, int, str, int]]:
        """
        Find tandem repeats.
        
        Args:
            sequence: DNA sequence
            min_repeat_length: Minimum length of repeat unit
            min_repeats: Minimum number of repeats
            
        Returns:
            List of (start, end, repeat_unit, count) tuples
        """
        sequence = sequence.upper()
        repeats = []
        
        for unit_len in range(min_repeat_length, len(sequence) // min_repeats + 1):
            for i in range(len(sequence) - unit_len * min_repeats + 1):
                unit = sequence[i:i + unit_len]
                count = 1
                
                # Count consecutive repeats
                j = i + unit_len
                while j + unit_len <= len(sequence):
                    if sequence[j:j + unit_len] == unit:
                        count += 1
                        j += unit_len
                    else:
                        break
                
                if count >= min_repeats:
                    end = i + count * unit_len
                    repeats.append((i, end, unit, count))
        
        # Remove overlapping repeats (keep longest)
        repeats.sort(key=lambda x: x[1] - x[0], reverse=True)
        filtered_repeats = []
        used_positions = set()
        
        for start, end, unit, count in repeats:
            positions = set(range(start, end))
            if not positions.intersection(used_positions):
                filtered_repeats.append((start, end, unit, count))
                used_positions.update(positions)
        
        return filtered_repeats
    
    def find_inverted_repeats(self,
                              sequence: str,
                              min_length: int = 10,
                              max_gap: int = 0) -> List[Tuple[int, int, int, int]]:
        """
        Find inverted repeats (palindromes).
        
        Args:
            sequence: DNA sequence
            min_length: Minimum length of repeat
            max_gap: Maximum gap between repeats
            
        Returns:
            List of (start1, end1, start2, end2) tuples
        """
        sequence = sequence.upper()
        repeats = []
        
        # Complement mapping
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        
        for i in range(len(sequence) - min_length + 1):
            for j in range(i + min_length, len(sequence) - min_length + 1):
                # Check if reverse complement matches
                seq1 = sequence[i:i + min_length]
                seq2 = sequence[j:j + min_length]
                rev_comp_seq2 = ''.join(complement.get(nuc, nuc) for nuc in reversed(seq2))
                
                if seq1 == rev_comp_seq2:
                    # Extend match
                    length = min_length
                    while (i + length < len(sequence) and 
                           j + length < len(sequence) and
                           sequence[i + length] == complement.get(sequence[j + length - 1], 'N')):
                        length += 1
                    
                    if length >= min_length:
                        gap = j - (i + length)
                        if gap <= max_gap:
                            repeats.append((i, i + length, j, j + length))
        
        return repeats
    
    def find_at_rich_regions(self,
                            sequence: str,
                            min_length: int = 20,
                            at_threshold: float = 0.7) -> List[Tuple[int, int, float]]:
        """
        Find AT-rich regions.
        
        Args:
            sequence: DNA sequence
            min_length: Minimum length of region
            at_threshold: Minimum AT content (0-1)
            
        Returns:
            List of (start, end, at_content) tuples
        """
        sequence = sequence.upper()
        regions = []
        
        for i in range(len(sequence) - min_length + 1):
            window = sequence[i:i + min_length]
            at_count = window.count('A') + window.count('T')
            at_content = at_count / len(window)
            
            if at_content >= at_threshold:
                # Extend region
                end = i + min_length
                while end < len(sequence):
                    if sequence[end] in 'AT':
                        end += 1
                    else:
                        break
                
                at_content_full = (sequence[i:end].count('A') + sequence[i:end].count('T')) / (end - i)
                regions.append((i, end, at_content_full))
        
        # Merge overlapping regions
        if not regions:
            return []
        
        regions.sort()
        merged = [regions[0]]
        
        for start, end, at_content in regions[1:]:
            last_start, last_end, _ = merged[-1]
            if start <= last_end:
                # Merge
                merged[-1] = (last_start, max(end, last_end), at_content)
            else:
                merged.append((start, end, at_content))
        
        return merged
    
    def find_gc_rich_regions(self,
                            sequence: str,
                            min_length: int = 20,
                            gc_threshold: float = 0.8) -> List[Tuple[int, int, float]]:
        """
        Find GC-rich regions.
        
        Args:
            sequence: DNA sequence
            min_length: Minimum length of region
            gc_threshold: Minimum GC content (0-1)
            
        Returns:
            List of (start, end, gc_content) tuples
        """
        sequence = sequence.upper()
        regions = []
        
        for i in range(len(sequence) - min_length + 1):
            window = sequence[i:i + min_length]
            gc_count = window.count('G') + window.count('C')
            gc_content = gc_count / len(window)
            
            if gc_content >= gc_threshold:
                # Extend region
                end = i + min_length
                while end < len(sequence):
                    if sequence[end] in 'GC':
                        end += 1
                    else:
                        break
                
                gc_content_full = (sequence[i:end].count('G') + sequence[i:end].count('C')) / (end - i)
                regions.append((i, end, gc_content_full))
        
        # Merge overlapping regions
        if not regions:
            return []
        
        regions.sort()
        merged = [regions[0]]
        
        for start, end, gc_content in regions[1:]:
            last_start, last_end, _ = merged[-1]
            if start <= last_end:
                merged[-1] = (last_start, max(end, last_end), gc_content)
            else:
                merged.append((start, end, gc_content))
        
        return merged
    
    def get_all_problematic_motifs(self, sequence: str) -> Dict[str, List]:
        """
        Find all problematic motifs in sequence.
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with all detected motifs
        """
        return {
            'tata_boxes': self.find_tata_boxes(sequence),
            'cryptic_promoters': self.find_cryptic_promoters(sequence),
            'splice_sites': self.find_splice_sites(sequence),
            'polyA_signals': self.find_polyA_signals(sequence),
            'homopolymers': self.find_homopolymers(sequence),
            'tandem_repeats': self.find_tandem_repeats(sequence),
            'inverted_repeats': self.find_inverted_repeats(sequence),
            'at_rich_regions': self.find_at_rich_regions(sequence),
            'gc_rich_regions': self.find_gc_rich_regions(sequence),
        }
    
    def calculate_motif_penalty(self, sequence: str) -> float:
        """
        Calculate penalty score for problematic motifs (0-1, lower is better).
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Penalty score
        """
        motifs = self.get_all_problematic_motifs(sequence)
        
        total_penalty = 0.0
        
        # Weight different motif types
        weights = {
            'tata_boxes': 0.2,
            'cryptic_promoters': 0.15,
            'splice_sites': 0.15,
            'polyA_signals': 0.15,
            'homopolymers': 0.1,
            'tandem_repeats': 0.1,
            'inverted_repeats': 0.1,
            'at_rich_regions': 0.05,
            'gc_rich_regions': 0.05,
        }
        
        for motif_type, weight in weights.items():
            motif_list = motifs.get(motif_type, [])
            if isinstance(motif_list, dict):
                # For splice_sites which is a dict
                count = sum(len(v) for v in motif_list.values())
            else:
                count = len(motif_list)
            
            # Normalize penalty by sequence length
            normalized_count = count / max(len(sequence) / 100, 1)
            total_penalty += weight * min(normalized_count, 1.0)
        
        return min(total_penalty, 1.0)




