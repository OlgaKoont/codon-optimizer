"""
CQA (Critical Quality Attributes) predictor and ranker for antibodies.
"""

from typing import List, Dict, Optional, Tuple
from codon_optimizer.core.sequence import OptimizedSequence
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from codon_optimizer.analysis.gc_content import GCContentAnalyzer
from codon_optimizer.constraints.protein_integrity import ProteinIntegrityChecker
from Bio.Seq import Seq


class CQAPredictor:
    """Predictor for Critical Quality Attributes of antibodies."""
    
    def __init__(self):
        """Initialize CQA predictor."""
        self.integrity_checker = ProteinIntegrityChecker()
        self.codon_analyzer = CodonUsageAnalyzer()
        self.gc_analyzer = GCContentAnalyzer()
    
    def predict_glycosylation_risk(self, sequence: OptimizedSequence) -> Dict[str, any]:
        """
        Predict glycosylation-related risks.
        
        Args:
            sequence: Optimized sequence
            
        Returns:
            Dictionary with glycosylation predictions
        """
        protein = sequence.get_protein_sequence()
        
        # Find N-glycosylation sites
        glycosylation_sites = self.integrity_checker.find_glycosylation_sites(protein)
        
        # Check for atypical sites in variable domains
        # (This is simplified - real analysis would identify CDR regions)
        atypical_sites = []
        for pos, motif in glycosylation_sites:
            # Check if in variable region (first ~120 amino acids for HC, ~110 for LC)
            if pos < 120:
                atypical_sites.append((pos, motif))
        
        return {
            'total_glycosylation_sites': len(glycosylation_sites),
            'atypical_sites': atypical_sites,
            'atypical_site_count': len(atypical_sites),
            'risk_level': 'high' if len(atypical_sites) > 0 else 'low'
        }
    
    def predict_aggregation_risk(self, sequence: OptimizedSequence) -> Dict[str, any]:
        """
        Predict aggregation risk based on sequence properties.
        
        Args:
            sequence: Optimized sequence
            
        Returns:
            Dictionary with aggregation predictions
        """
        protein = sequence.get_protein_sequence()
        
        # Find hydrophobic regions
        hydrophobic_regions = self.integrity_checker.find_hydrophobic_regions(
            protein, min_length=10
        )
        
        # Find unpaired cysteines
        unpaired_cys = self.integrity_checker.find_unpaired_cysteines(protein)
        
        # Calculate hydrophobicity index (simplified)
        hydrophobic_aas = set('AILMVFYW')
        hydrophobic_count = sum(1 for aa in protein if aa in hydrophobic_aas)
        hydrophobicity_index = hydrophobic_count / len(protein) if protein else 0
        
        # Risk assessment
        risk_factors = []
        if len(hydrophobic_regions) > 2:
            risk_factors.append('multiple_hydrophobic_regions')
        if len(unpaired_cys) > 0:
            risk_factors.append('unpaired_cysteines')
        if hydrophobicity_index > 0.4:
            risk_factors.append('high_hydrophobicity')
        
        risk_level = 'high' if len(risk_factors) >= 2 else 'medium' if len(risk_factors) == 1 else 'low'
        
        return {
            'hydrophobic_regions': len(hydrophobic_regions),
            'unpaired_cysteines': len(unpaired_cys),
            'hydrophobicity_index': hydrophobicity_index,
            'risk_factors': risk_factors,
            'risk_level': risk_level
        }
    
    def predict_charge_variant_risk(self, sequence: OptimizedSequence) -> Dict[str, any]:
        """
        Predict charge variant risk.
        
        Args:
            sequence: Optimized sequence
            
        Returns:
            Dictionary with charge variant predictions
        """
        protein = sequence.get_protein_sequence()
        
        # Count charged amino acids
        positive_charged = set('KRH')
        negative_charged = set('DE')
        
        positive_count = sum(1 for aa in protein if aa in positive_charged)
        negative_count = sum(1 for aa in protein if aa in negative_charged)
        
        net_charge = positive_count - negative_count
        charge_density = (positive_count + negative_count) / len(protein) if protein else 0
        
        # Risk: high charge density or extreme net charge
        risk_level = 'low'
        if abs(net_charge) > 20:
            risk_level = 'high'
        elif abs(net_charge) > 10:
            risk_level = 'medium'
        
        return {
            'positive_charged': positive_count,
            'negative_charged': negative_count,
            'net_charge': net_charge,
            'charge_density': charge_density,
            'risk_level': risk_level
        }
    
    def predict_expression_level(self, sequence: OptimizedSequence) -> Dict[str, any]:
        """
        Predict expression level based on codon usage.
        
        Args:
            sequence: Optimized sequence
            
        Returns:
            Dictionary with expression predictions
        """
        cai = self.codon_analyzer.calculate_cai(sequence.dna_sequence)
        
        # Estimate expression level (simplified model)
        # In practice, this would use experimental data
        if cai >= 0.8:
            expression_level = 'high'
        elif cai >= 0.6:
            expression_level = 'medium'
        else:
            expression_level = 'low'
        
        # Check for tRNA depletion
        has_depletion, usage_ratios = self.codon_analyzer.check_tRNA_depletion(
            sequence.dna_sequence
        )
        
        if has_depletion:
            expression_level = 'reduced'
        
        return {
            'cai': cai,
            'predicted_expression': expression_level,
            'tRNA_depletion_risk': has_depletion,
            'max_codon_usage': max(usage_ratios.values()) if usage_ratios else 0.0
        }
    
    def get_all_cqa_predictions(self, sequence: OptimizedSequence) -> Dict[str, any]:
        """
        Get all CQA predictions for a sequence.
        
        Args:
            sequence: Optimized sequence
            
        Returns:
            Dictionary with all CQA predictions
        """
        return {
            'glycosylation': self.predict_glycosylation_risk(sequence),
            'aggregation': self.predict_aggregation_risk(sequence),
            'charge_variants': self.predict_charge_variant_risk(sequence),
            'expression': self.predict_expression_level(sequence)
        }
    
    def calculate_cqa_score(self, sequence: OptimizedSequence) -> float:
        """
        Calculate overall CQA score (0-1, higher is better).
        
        Args:
            sequence: Optimized sequence
            
        Returns:
            CQA score
        """
        predictions = self.get_all_cqa_predictions(sequence)
        
        score = 1.0
        
        # Penalties for risks
        if predictions['glycosylation']['risk_level'] == 'high':
            score -= 0.2
        elif predictions['glycosylation']['risk_level'] == 'medium':
            score -= 0.1
        
        if predictions['aggregation']['risk_level'] == 'high':
            score -= 0.3
        elif predictions['aggregation']['risk_level'] == 'medium':
            score -= 0.15
        
        if predictions['charge_variants']['risk_level'] == 'high':
            score -= 0.2
        elif predictions['charge_variants']['risk_level'] == 'medium':
            score -= 0.1
        
        if predictions['expression']['predicted_expression'] == 'low':
            score -= 0.2
        elif predictions['expression']['predicted_expression'] == 'reduced':
            score -= 0.15
        
        return max(0.0, score)
    
    def rank_sequences(self,
                      sequences: List[OptimizedSequence],
                      experimental_cqa: Optional[Dict[str, Dict]] = None) -> List[Tuple[OptimizedSequence, float, Dict]]:
        """
        Rank sequences based on CQA predictions and optionally experimental data.
        
        Args:
            sequences: List of sequences to rank
            experimental_cqa: Optional dictionary mapping sequence IDs to experimental CQA data
            
        Returns:
            List of (sequence, cqa_score, predictions) tuples, sorted by score
        """
        ranked = []
        
        for seq in sequences:
            predictions = self.get_all_cqa_predictions(seq)
            cqa_score = self.calculate_cqa_score(seq)
            
            # Adjust score based on experimental data if available
            if experimental_cqa and seq.sequence_id in experimental_cqa:
                exp_data = experimental_cqa[seq.sequence_id]
                
                # Weight experimental data more heavily
                exp_score = self._calculate_experimental_score(exp_data)
                cqa_score = 0.7 * cqa_score + 0.3 * exp_score
            
            ranked.append((seq, cqa_score, predictions))
        
        # Sort by score (descending)
        ranked.sort(key=lambda x: x[1], reverse=True)
        
        return ranked
    
    def _calculate_experimental_score(self, exp_data: Dict) -> float:
        """
        Calculate score from experimental CQA data.
        
        Args:
            exp_data: Experimental CQA data dictionary
            
        Returns:
            Score (0-1)
        """
        score = 1.0
        
        # Penalize based on experimental issues
        if exp_data.get('aggregation', {}).get('high', False):
            score -= 0.3
        if exp_data.get('glycosylation', {}).get('atypical', False):
            score -= 0.2
        if exp_data.get('charge_variants', {}).get('high', False):
            score -= 0.2
        if exp_data.get('expression', {}).get('low', False):
            score -= 0.2
        
        return max(0.0, score)
    
    def filter_by_cqa(self,
                      sequences: List[OptimizedSequence],
                      min_cqa_score: float = 0.6,
                      exclude_high_risk: bool = True) -> List[OptimizedSequence]:
        """
        Filter sequences based on CQA criteria.
        
        Args:
            sequences: List of sequences to filter
            min_cqa_score: Minimum CQA score
            exclude_high_risk: Exclude sequences with high-risk CQA issues
            
        Returns:
            Filtered list of sequences
        """
        filtered = []
        
        for seq in sequences:
            predictions = self.get_all_cqa_predictions(seq)
            cqa_score = self.calculate_cqa_score(seq)
            
            # Check minimum score
            if cqa_score < min_cqa_score:
                continue
            
            # Check for high-risk issues
            if exclude_high_risk:
                has_high_risk = (
                    predictions['glycosylation']['risk_level'] == 'high' or
                    predictions['aggregation']['risk_level'] == 'high' or
                    predictions['charge_variants']['risk_level'] == 'high'
                )
                if has_high_risk:
                    continue
            
            filtered.append(seq)
        
        return filtered




