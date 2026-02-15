"""
Report generator for codon optimization results.
"""

import json
from typing import List, Dict, Optional
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from codon_optimizer.core.sequence import OptimizedSequence, SequencePair
from codon_optimizer.core.scorer import MultiCriteriaScorer
from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
from codon_optimizer.analysis.gc_content import GCContentAnalyzer
from codon_optimizer.antibody.cqa_predictor import CQAPredictor


class ReportGenerator:
    """Generator for optimization reports."""
    
    def __init__(self, scorer: MultiCriteriaScorer):
        """
        Initialize report generator.
        
        Args:
            scorer: MultiCriteriaScorer instance
        """
        self.scorer = scorer
        self.codon_analyzer = CodonUsageAnalyzer()
        self.gc_analyzer = GCContentAnalyzer()
        self.cqa_predictor = CQAPredictor()
    
    def generate_text_report(self,
                            sequences: List[OptimizedSequence],
                            output_path: str,
                            original_sequence: Optional[str] = None):
        """
        Generate text report.
        
        Args:
            sequences: List of optimized sequences
            original_sequence: Original DNA sequence
            output_path: Path to output file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CODON OPTIMIZATION REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            if original_sequence:
                f.write(f"Original Sequence Length: {len(original_sequence)} bp\n")
                f.write(f"Original Protein Length: {len(str(Seq(original_sequence).translate()))} aa\n\n")
            
            f.write(f"Number of Optimized Solutions: {len(sequences)}\n\n")
            
            for i, seq in enumerate(sequences, 1):
                f.write("-" * 80 + "\n")
                f.write(f"Solution {i}: {seq.sequence_id}\n")
                f.write("-" * 80 + "\n\n")
                
                # Get detailed analysis
                original_protein = seq.get_protein_sequence() if seq.original_sequence else None
                analysis = self.scorer.get_detailed_analysis(seq.dna_sequence, original_protein)
                
                # Scores
                scores = analysis['scores']
                f.write("SCORES:\n")
                f.write(f"  Total Score: {scores['total_score']:.4f}\n")
                f.write(f"  CAI Score: {scores['cai_score']:.4f}\n")
                f.write(f"  GC Score: {scores['gc_score']:.4f}\n")
                f.write(f"  Pair Score: {scores['pair_score']:.4f}\n")
                f.write(f"  Structure Score: {scores['structure_score']:.4f}\n")
                f.write(f"  Motif Score: {scores['motif_score']:.4f}\n")
                f.write(f"  Cloning Score: {scores['cloning_score']:.4f}\n\n")
                
                # CAI
                f.write(f"CAI: {analysis['cai']:.4f}\n\n")
                
                # GC Statistics
                gc_stats = analysis['gc_statistics']
                f.write("GC CONTENT:\n")
                f.write(f"  Overall GC: {gc_stats['overall_gc']:.4f}\n")
                f.write(f"  Mean GC: {gc_stats['mean_gc']:.4f}\n")
                f.write(f"  Min GC: {gc_stats['min_gc']:.4f}\n")
                f.write(f"  Max GC: {gc_stats['max_gc']:.4f}\n")
                f.write(f"  Violations: {gc_stats['violation_count']}\n\n")
                
                # Structure
                struct_stats = analysis['structure_statistics']
                f.write("mRNA STRUCTURE:\n")
                f.write(f"  5' Region ΔG: {struct_stats.get('critical_region_dg', 'N/A')}\n")
                f.write(f"  Structure Score: {struct_stats.get('structure_score', 0.0):.4f}\n")
                f.write(f"  Stable Structures: {struct_stats.get('stable_structure_count', 0)}\n\n")
                
                # Motifs
                motifs = analysis['motifs']
                f.write("PROBLEMATIC MOTIFS:\n")
                for motif_type, motif_list in motifs.items():
                    if isinstance(motif_list, dict):
                        count = sum(len(v) for v in motif_list.values())
                    else:
                        count = len(motif_list)
                    if count > 0:
                        f.write(f"  {motif_type}: {count}\n")
                f.write("\n")
                
                # Restriction Sites
                restriction = analysis['restriction_sites']
                f.write("RESTRICTION SITES:\n")
                f.write(f"  Unwanted Sites: {restriction.get('violation_count', 0)}\n")
                if restriction.get('missing_required_sites'):
                    f.write(f"  Missing Required: {', '.join(restriction['missing_required_sites'])}\n")
                f.write("\n")
                
                # CQA Predictions
                cqa_predictions = self.cqa_predictor.get_all_cqa_predictions(seq)
                f.write("CQA PREDICTIONS:\n")
                f.write(f"  Glycosylation Risk: {cqa_predictions['glycosylation']['risk_level']}\n")
                f.write(f"  Aggregation Risk: {cqa_predictions['aggregation']['risk_level']}\n")
                f.write(f"  Charge Variant Risk: {cqa_predictions['charge_variants']['risk_level']}\n")
                f.write(f"  Expression Level: {cqa_predictions['expression']['predicted_expression']}\n\n")
                
                f.write("\n")
    
    def generate_json_report(self,
                            sequences: List[OptimizedSequence],
                            output_path: str):
        """
        Generate JSON report.
        
        Args:
            sequences: List of optimized sequences
            output_path: Path to output file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        report_data = {
            'solutions': []
        }
        
        for seq in sequences:
            original_protein = seq.get_protein_sequence() if seq.original_sequence else None
            analysis = self.scorer.get_detailed_analysis(seq.dna_sequence, original_protein)
            cqa_predictions = self.cqa_predictor.get_all_cqa_predictions(seq)
            
            solution_data = {
                'sequence_id': seq.sequence_id,
                'dna_sequence': seq.dna_sequence,
                'protein_sequence': seq.get_protein_sequence(),
                'scores': analysis['scores'],
                'cai': analysis['cai'],
                'gc_statistics': analysis['gc_statistics'],
                'pair_statistics': analysis['pair_statistics'],
                'structure_statistics': analysis['structure_statistics'],
                'motifs': {k: len(v) if isinstance(v, list) else v for k, v in analysis['motifs'].items()},
                'restriction_sites': analysis['restriction_sites'],
                'cqa_predictions': cqa_predictions,
                'metadata': seq.metadata
            }
            
            report_data['solutions'].append(solution_data)
        
        with open(output_path, 'w') as f:
            json.dump(report_data, f, indent=2)
    
    def generate_html_report(self,
                            sequences: List[OptimizedSequence],
                            output_path: str,
                            original_sequence: Optional[str] = None):
        """
        Generate HTML report with visualizations.
        
        Args:
            sequences: List of optimized sequences
            output_path: Path to output file
            original_sequence: Original DNA sequence
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        html = []
        html.append("<!DOCTYPE html>")
        html.append("<html><head>")
        html.append("<title>Codon Optimization Report</title>")
        html.append("<style>")
        html.append("body { font-family: Arial, sans-serif; margin: 20px; }")
        html.append("h1 { color: #333; }")
        html.append("h2 { color: #666; border-bottom: 2px solid #ccc; padding-bottom: 5px; }")
        html.append("table { border-collapse: collapse; width: 100%; margin: 20px 0; }")
        html.append("th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }")
        html.append("th { background-color: #4CAF50; color: white; }")
        html.append("tr:nth-child(even) { background-color: #f2f2f2; }")
        html.append(".score { font-weight: bold; }")
        html.append(".high { color: green; }")
        html.append(".medium { color: orange; }")
        html.append(".low { color: red; }")
        html.append("</style>")
        html.append("</head><body>")
        
        html.append("<h1>Codon Optimization Report</h1>")
        
        if original_sequence:
            html.append(f"<p><strong>Original Sequence Length:</strong> {len(original_sequence)} bp</p>")
        
        html.append(f"<p><strong>Number of Solutions:</strong> {len(sequences)}</p>")
        
        # Summary table
        html.append("<h2>Summary</h2>")
        html.append("<table>")
        html.append("<tr><th>Solution</th><th>Total Score</th><th>CAI</th><th>GC Score</th><th>Structure Score</th><th>CQA Score</th></tr>")
        
        for seq in sequences:
            original_protein = seq.get_protein_sequence() if seq.original_sequence else None
            analysis = self.scorer.get_detailed_analysis(seq.dna_sequence, original_protein)
            cqa_score = self.cqa_predictor.calculate_cqa_score(seq)
            
            scores = analysis['scores']
            html.append(f"<tr>")
            html.append(f"<td>{seq.sequence_id}</td>")
            html.append(f"<td class='score'>{scores['total_score']:.4f}</td>")
            html.append(f"<td>{analysis['cai']:.4f}</td>")
            html.append(f"<td>{scores['gc_score']:.4f}</td>")
            html.append(f"<td>{scores['structure_score']:.4f}</td>")
            html.append(f"<td>{cqa_score:.4f}</td>")
            html.append(f"</tr>")
        
        html.append("</table>")
        
        # Detailed sections for each solution
        for i, seq in enumerate(sequences, 1):
            html.append(f"<h2>Solution {i}: {seq.sequence_id}</h2>")
            
            original_protein = seq.get_protein_sequence() if seq.original_sequence else None
            analysis = self.scorer.get_detailed_analysis(seq.dna_sequence, original_protein)
            cqa_predictions = self.cqa_predictor.get_all_cqa_predictions(seq)
            
            # Scores
            scores = analysis['scores']
            html.append("<h3>Scores</h3>")
            html.append("<table>")
            html.append("<tr><th>Metric</th><th>Value</th></tr>")
            html.append(f"<tr><td>Total Score</td><td class='score'>{scores['total_score']:.4f}</td></tr>")
            html.append(f"<tr><td>CAI Score</td><td>{scores['cai_score']:.4f}</td></tr>")
            html.append(f"<tr><td>GC Score</td><td>{scores['gc_score']:.4f}</td></tr>")
            html.append(f"<tr><td>Pair Score</td><td>{scores['pair_score']:.4f}</td></tr>")
            html.append(f"<tr><td>Structure Score</td><td>{scores['structure_score']:.4f}</td></tr>")
            html.append(f"<tr><td>Motif Score</td><td>{scores['motif_score']:.4f}</td></tr>")
            html.append(f"<tr><td>Cloning Score</td><td>{scores['cloning_score']:.4f}</td></tr>")
            html.append("</table>")
            
            # CQA Predictions
            html.append("<h3>CQA Predictions</h3>")
            html.append("<table>")
            html.append("<tr><th>Attribute</th><th>Risk Level</th><th>Details</th></tr>")
            
            for attr, pred in cqa_predictions.items():
                risk_level = pred.get('risk_level', pred.get('predicted_expression', 'N/A'))
                risk_class = 'high' if risk_level in ['high', 'low', 'reduced'] else 'medium' if risk_level == 'medium' else 'low'
                html.append(f"<tr>")
                html.append(f"<td>{attr.capitalize()}</td>")
                html.append(f"<td class='{risk_class}'>{risk_level}</td>")
                html.append(f"<td>{json.dumps(pred, indent=2)}</td>")
                html.append(f"</tr>")
            
            html.append("</table>")
            
            # Sequence
            html.append("<h3>Optimized Sequence</h3>")
            html.append(f"<p><code>{seq.dna_sequence}</code></p>")
        
        html.append("</body></html>")
        
        with open(output_path, 'w') as f:
            f.write('\n'.join(html))
    
    def plot_gc_profile(self, sequence: OptimizedSequence, output_path: str):
        """
        Plot GC content profile.
        
        Args:
            sequence: Optimized sequence
            output_path: Path to output image
        """
        gc_profile = self.gc_analyzer.calculate_gc_profile(sequence.dna_sequence)
        
        plt.figure(figsize=(10, 6))
        plt.plot(gc_profile, label='GC Content')
        plt.axhline(y=self.gc_analyzer.gc_min, color='r', linestyle='--', label='Min GC')
        plt.axhline(y=self.gc_analyzer.gc_max, color='r', linestyle='--', label='Max GC')
        plt.xlabel('Position (window)')
        plt.ylabel('GC Content')
        plt.title(f'GC Content Profile - {sequence.sequence_id}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()
    
    def plot_comparison(self,
                       sequences: List[OptimizedSequence],
                       output_path: str):
        """
        Plot comparison of multiple solutions.
        
        Args:
            sequences: List of sequences to compare
            output_path: Path to output image
        """
        metrics = ['cai', 'gc_score', 'structure_score', 'motif_score']
        data = {metric: [] for metric in metrics}
        labels = []
        
        for seq in sequences:
            labels.append(seq.sequence_id)
            original_protein = seq.get_protein_sequence() if seq.original_sequence else None
            analysis = self.scorer.get_detailed_analysis(seq.dna_sequence, original_protein)
            
            data['cai'].append(analysis['cai'])
            scores = analysis['scores']
            data['gc_score'].append(scores['gc_score'])
            data['structure_score'].append(scores['structure_score'])
            data['motif_score'].append(scores['motif_score'])
        
        x = np.arange(len(labels))
        width = 0.2
        
        fig, ax = plt.subplots(figsize=(12, 6))
        for i, metric in enumerate(metrics):
            offset = (i - len(metrics)/2) * width
            ax.bar(x + offset, data[metric], width, label=metric.replace('_', ' ').title())
        
        ax.set_xlabel('Solution')
        ax.set_ylabel('Score')
        ax.set_title('Comparison of Optimized Solutions')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()




