"""
Command-line interface for codon optimization.
"""

import click
import sys
import traceback
from pathlib import Path
from codon_optimizer.core.config import OptimizationConfig
from codon_optimizer.core.optimizer import CodonOptimizer
from codon_optimizer.core.pareto import ParetoOptimizer
from codon_optimizer.core.scorer import MultiCriteriaScorer
from codon_optimizer.io.fasta_handler import FastaHandler
from codon_optimizer.io.report_generator import ReportGenerator
from codon_optimizer.antibody.antibody_handler import AntibodyOptimizer
from codon_optimizer.antibody.cqa_predictor import CQAPredictor


@click.group()
def cli():
    """Codon Optimization System for Antibody Production."""
    pass


@cli.command()
@click.option('--input', '-i', required=True, help='Input FASTA file')
@click.option('--output', '-o', required=True, help='Output FASTA file')
@click.option('--host', default='CHO', help='Host organism (default: CHO)')
@click.option('--mode', type=click.Choice(['standard', 'pareto', 'antibody']), 
              default='standard', help='Optimization mode')
@click.option('--heavy-chain', help='Heavy chain FASTA file (for antibody mode)')
@click.option('--light-chain', help='Light chain FASTA file (for antibody mode)')
@click.option('--restriction-sites', help='Comma-separated list of restriction sites to remove')
@click.option('--target-expression', help='Target expression range (e.g., 0.7-0.9)')
@click.option('--pareto-solutions', default=10, help='Number of Pareto solutions (default: 10)')
@click.option('--num-solutions', default=1, help='Number of diverse solutions to generate in standard mode (default: 1, max: 5)')
@click.option('--report', help='Generate report (text, json, or html)')
@click.option('--report-format', type=click.Choice(['text', 'json', 'html']), 
              default='html', help='Report format (default: html)')
@click.option('--population-size', default=50, help='Population size for GA (default: 50)')
@click.option('--max-generations', default=100, help='Max generations for GA (default: 100)')
def optimize(input, output, host, mode, heavy_chain, light_chain, restriction_sites,
            target_expression, pareto_solutions, num_solutions, report, report_format, 
            population_size, max_generations):
    """Optimize codon usage for a DNA sequence."""
    
    # Parse restriction sites
    restriction_sites_list = []
    if restriction_sites:
        restriction_sites_list = [s.strip() for s in restriction_sites.split(',')]
    
    # Parse target expression
    target_expr_min, target_expr_max = 0.7, 0.9
    if target_expression:
        parts = target_expression.split('-')
        if len(parts) == 2:
            target_expr_min = float(parts[0])
            target_expr_max = float(parts[1])
    
    # Create configuration
    config = OptimizationConfig(
        host=host,
        restriction_sites_to_remove=restriction_sites_list,
        target_expression_min=target_expr_min,
        target_expression_max=target_expr_max,
        population_size=population_size,
        max_generations=max_generations,
        use_pareto=(mode == 'pareto'),
        pareto_solutions=pareto_solutions,
        is_antibody=(mode == 'antibody')
    )
    try:
        config.validate()
    except ValueError as e:
        click.echo(f"Configuration error: {e}", err=True)
        sys.exit(1)
    
    # Read input sequence - detect if it's protein or DNA
    try:
        seq_id, sequence = FastaHandler.read_single_sequence(input)
        sequence_type = FastaHandler.detect_sequence_type(sequence)
        
        # Convert DNA to protein if needed
        if sequence_type == "DNA":
            click.echo("Input detected as DNA sequence. Converting to protein...")
            protein_sequence = FastaHandler.translate_dna_to_protein(sequence)
        elif sequence_type == "PROTEIN":
            click.echo("Input detected as protein sequence.")
            protein_sequence = sequence
        else:
            click.echo(f"Error: Unsupported sequence type: {sequence_type}", err=True)
            sys.exit(1)
            
    except KeyboardInterrupt:
        click.echo("\nOperation cancelled by user.", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"Error reading input file: {e}", err=True)
        click.echo(traceback.format_exc(), err=True)
        sys.exit(1)
    
    # Optimize based on mode
    if mode == 'antibody':
        if not heavy_chain or not light_chain:
            click.echo("Error: --heavy-chain and --light-chain required for antibody mode", err=True)
            sys.exit(1)
        
        try:
            hc_id, hc_seq = FastaHandler.read_single_sequence(heavy_chain)
            lc_id, lc_seq = FastaHandler.read_single_sequence(light_chain)
            
            # Detect and convert to protein if needed
            hc_type = FastaHandler.detect_sequence_type(hc_seq)
            lc_type = FastaHandler.detect_sequence_type(lc_seq)
            
            if hc_type == "DNA":
                hc_protein = FastaHandler.translate_dna_to_protein(hc_seq)
            else:
                hc_protein = hc_seq
                
            if lc_type == "DNA":
                lc_protein = FastaHandler.translate_dna_to_protein(lc_seq)
            else:
                lc_protein = lc_seq
                
        except KeyboardInterrupt:
            click.echo("\nOperation cancelled by user.", err=True)
            sys.exit(130)
        except Exception as e:
            click.echo(f"Error reading antibody chains: {e}", err=True)
            click.echo(traceback.format_exc(), err=True)
            sys.exit(1)
        
        try:
            antibody_optimizer = AntibodyOptimizer(config)
            # Note: antibody optimizer needs to be updated to work with proteins
            # For now, we'll optimize each chain separately
            optimizer = CodonOptimizer(config)
            hc_optimized = optimizer.optimize(hc_protein)
            lc_optimized = optimizer.optimize(lc_protein)
            
            from codon_optimizer.core.sequence import SequencePair
            pair = SequencePair(hc_optimized, lc_optimized, "HC_LC")
            
            # Write outputs
            FastaHandler.write_fasta([
                (f"{hc_id}_optimized", pair.sequence1.dna_sequence),
                (f"{lc_id}_optimized", pair.sequence2.dna_sequence)
            ], output)
            
            click.echo(f"Optimized antibody chains written to {output}")
            
            # Generate report if requested
            if report:
                scorer = MultiCriteriaScorer(config)
                report_gen = ReportGenerator(scorer)
                
                if report_format == 'html':
                    report_gen.generate_html_report(
                        [pair.sequence1, pair.sequence2],
                        report,
                        original_sequence=None
                    )
                elif report_format == 'json':
                    report_gen.generate_json_report([pair.sequence1, pair.sequence2], report)
                else:
                    report_gen.generate_text_report(
                        [pair.sequence1, pair.sequence2],
                        report,
                        original_sequence=None
                    )
                
                click.echo(f"Report written to {report}")
        except KeyboardInterrupt:
            click.echo("\nOperation cancelled by user during antibody optimization.", err=True)
            sys.exit(130)
        except Exception as e:
            error_msg = f"Error during antibody optimization: {str(e)}"
            click.echo(error_msg, err=True)
            click.echo("\nFull traceback:", err=True)
            click.echo(traceback.format_exc(), err=True)
            sys.exit(1)
    
    elif mode == 'pareto':
        try:
            # For Pareto, use ParetoOptimizer
            pareto_optimizer = ParetoOptimizer(config)
            solutions = pareto_optimizer.optimize_pareto(protein_sequence, pareto_solutions)
            
            # Write all solutions
            output_sequences = [
                (f"{seq_id}_optimized_{i+1}", sol.dna_sequence)
                for i, sol in enumerate(solutions)
            ]
            FastaHandler.write_fasta(output_sequences, output)
            
            click.echo(f"Generated {len(solutions)} optimized solutions")
            click.echo(f"Solutions written to {output}")
            
            # Generate report
            if report:
                scorer = MultiCriteriaScorer(config)
                report_gen = ReportGenerator(scorer)
                
                if report_format == 'html':
                    report_gen.generate_html_report(solutions, report, original_sequence=None)
                elif report_format == 'json':
                    report_gen.generate_json_report(solutions, report)
                else:
                    report_gen.generate_text_report(solutions, report, original_sequence=None)
                
                click.echo(f"Report written to {report}")
        except KeyboardInterrupt:
            click.echo("\nOperation cancelled by user during Pareto optimization.", err=True)
            sys.exit(130)
        except Exception as e:
            error_msg = f"Error during Pareto optimization: {str(e)}"
            click.echo(error_msg, err=True)
            click.echo("\nFull traceback:", err=True)
            click.echo(traceback.format_exc(), err=True)
            sys.exit(1)
    
    else:  # standard mode
        try:
            optimizer = CodonOptimizer(config)
            click.echo(f"Optimizing protein sequence (length: {len(protein_sequence)} aa)...")
            
            # Progress callback for real-time updates
            def progress_log(msg):
                # Use err=True to ensure immediate output (stderr is unbuffered)
                click.echo(msg, err=True)
                import sys
                sys.stderr.flush()  # Force immediate flush
            
            # Generate single or multiple solutions
            if num_solutions > 1:
                click.echo(f"Generating {num_solutions} diverse solutions (optimized for different criteria)...")
                # Use V2 optimizer for multiple diverse solutions
                from codon_optimizer.core.optimizer_v2 import CodonOptimizerV2
                optimizer_v2 = CodonOptimizerV2(config)
                solutions = optimizer_v2.optimize_multiple_diverse(
                    protein_sequence,
                    num_solutions=min(num_solutions, 5),  # Max 5 strategies
                    progress_callback=progress_log
                )
                
                # Write all solutions to FASTA
                output_sequences = [
                    (sol.sequence_id or f"{seq_id}_optimized_{i+1}", sol.dna_sequence)
                    for i, sol in enumerate(solutions)
                ]
                FastaHandler.write_fasta(output_sequences, output)
                
                click.echo(f"\nGenerated {len(solutions)} diverse solutions:")
                for i, sol in enumerate(solutions, 1):
                    strategy_desc = sol.metadata.get('strategy_description', 'optimized')
                    score = sol.metadata.get('fitness_score', 0)
                    cai = sol.metadata.get('cai', 0)
                    gc = sol.metadata.get('gc_score', 0)
                    click.echo(f"  {i}. {strategy_desc}: score={score:.4f}, CAI={cai:.4f}, GC={gc:.4f}")
                
                click.echo(f"\nAll solutions written to {output}")
                
                # Generate report for all solutions
                if report:
                    scorer = MultiCriteriaScorer(config)
                    report_gen = ReportGenerator(scorer)
                    
                    if report_format == 'html':
                        report_gen.generate_html_report(solutions, report, original_sequence=None)
                    elif report_format == 'json':
                        report_gen.generate_json_report(solutions, report)
                    else:
                        report_gen.generate_text_report(solutions, report, original_sequence=None)
                    
                    click.echo(f"Report written to {report}")
            else:
                # Single solution
                optimized = optimizer.optimize(protein_sequence, progress_callback=progress_log)
                
                # Show optimization method used
                method = optimized.metadata.get('optimization_method', 'unknown')
                click.echo(f"Optimization method: {method}")
                
                if 'total_combinations_evaluated' in optimized.metadata:
                    click.echo(f"Evaluated {optimized.metadata.get('total_combinations_evaluated', 0)} combinations")
                    click.echo(f"Total possible combinations: {optimized.metadata.get('total_possible_combinations', 0)}")
                
                if 'constraints_passed' in optimized.metadata:
                    click.echo(f"Constraints passed: {optimized.metadata.get('constraints_passed', False)}")
                    click.echo(f"Objectives score: {optimized.metadata.get('objectives_score', 0):.4f}")
                
                # Write output
                FastaHandler.write_single_sequence(
                    f"{seq_id}_optimized",
                    optimized.dna_sequence,
                    output
                )
                
                click.echo(f"Optimized sequence written to {output}")
                
                # Generate report
                if report:
                    scorer = MultiCriteriaScorer(config)
                    report_gen = ReportGenerator(scorer)
                    
                    if report_format == 'html':
                        report_gen.generate_html_report([optimized], report, original_sequence=None)
                    elif report_format == 'json':
                        report_gen.generate_json_report([optimized], report)
                    else:
                        report_gen.generate_text_report([optimized], report, original_sequence=None)
                    
                    click.echo(f"Report written to {report}")
        except KeyboardInterrupt:
            click.echo("\nOperation cancelled by user during optimization.", err=True)
            sys.exit(130)
        except Exception as e:
            error_msg = f"Error during optimization: {str(e)}"
            click.echo(error_msg, err=True)
            click.echo("\nFull traceback:", err=True)
            click.echo(traceback.format_exc(), err=True)
            sys.exit(1)


@cli.command()
@click.option('--input', '-i', required=True, help='Input FASTA file (protein or DNA sequence)')
@click.option('--output', '-o', help='Output file for analysis')
def analyze(input, output):
    """Analyze a sequence without optimization."""
    
    try:
        seq_id, sequence = FastaHandler.read_single_sequence(input)
        sequence_type = FastaHandler.detect_sequence_type(sequence)
        
        # Convert to protein if DNA
        if sequence_type == "DNA":
            click.echo("Input detected as DNA sequence. Converting to protein for analysis...")
            protein_sequence = FastaHandler.translate_dna_to_protein(sequence)
            # For analysis, we need a DNA sequence - generate one optimal
            from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
            codon_analyzer = CodonUsageAnalyzer()
            dna_sequence = ''.join([codon_analyzer.get_optimal_codon(aa) for aa in protein_sequence])
        elif sequence_type == "PROTEIN":
            click.echo("Input detected as protein sequence.")
            protein_sequence = sequence
            # Generate optimal DNA for analysis
            from codon_optimizer.analysis.codon_usage import CodonUsageAnalyzer
            codon_analyzer = CodonUsageAnalyzer()
            dna_sequence = ''.join([codon_analyzer.get_optimal_codon(aa) for aa in protein_sequence])
        else:
            click.echo(f"Error: Unsupported sequence type: {sequence_type}", err=True)
            sys.exit(1)
            
    except KeyboardInterrupt:
        click.echo("\nOperation cancelled by user.", err=True)
        sys.exit(130)
    except Exception as e:
        click.echo(f"Error reading input file: {e}", err=True)
        click.echo(traceback.format_exc(), err=True)
        sys.exit(1)
    
    try:
        config = OptimizationConfig()
        scorer = MultiCriteriaScorer(config)
        
        # Get detailed analysis
        analysis = scorer.get_detailed_analysis(dna_sequence, protein_sequence)
        
        # Print summary
        click.echo(f"Analysis for sequence: {seq_id}")
        click.echo(f"Sequence type: {sequence_type}")
        if sequence_type == "PROTEIN":
            click.echo(f"Protein length: {len(protein_sequence)} aa")
            click.echo(f"DNA length (optimal): {len(dna_sequence)} bp")
        else:
            click.echo(f"DNA length: {len(dna_sequence)} bp")
            click.echo(f"Protein length: {len(protein_sequence)} aa")
        click.echo(f"CAI: {analysis['cai']:.4f}")
        click.echo(f"Total Score: {analysis['scores']['total_score']:.4f}")
        click.echo(f"GC Content: {analysis['gc_statistics']['overall_gc']:.4f}")
        
        # Write detailed analysis if output specified
        if output:
            report_gen = ReportGenerator(scorer)
            from codon_optimizer.core.sequence import OptimizedSequence
            opt_seq = OptimizedSequence(dna_sequence, sequence_id=seq_id)
            
            if output.endswith('.json'):
                report_gen.generate_json_report([opt_seq], output)
            elif output.endswith('.html'):
                report_gen.generate_html_report([opt_seq], output, original_sequence=None)
            else:
                report_gen.generate_text_report([opt_seq], output, original_sequence=None)
            
            click.echo(f"Detailed analysis written to {output}")
    except KeyboardInterrupt:
        click.echo("\nOperation cancelled by user during analysis.", err=True)
        sys.exit(130)
    except Exception as e:
        error_msg = f"Error during analysis: {str(e)}"
        click.echo(error_msg, err=True)
        click.echo("\nFull traceback:", err=True)
        click.echo(traceback.format_exc(), err=True)
        sys.exit(1)


if __name__ == '__main__':
    cli()



