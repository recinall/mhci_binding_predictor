import os
import click
import pandas as pd
from .core import IEDBBindingPredictor

@click.group(epilog="""
EXAMPLES:
  Basic prediction:
  mhci_predict predict --peptides PEPTLIST.txt --alleles HLA-A*02:01,HLA-B*07:02
  
  Full analysis with pattern:
  mhci_predict full-analysis --pattern "A[CD]E[FY]GH" --alleles HLA-A*02:01 --output results.csv
  
  Filter results:
  mhci_predict filter --input results.csv --score 0.5 --ic50 500
""")
@click.option('--output-dir', default='./output', 
             help='Output directory for results (default: ./output)')
@click.option('--csv-sep', default=';', 
             help='CSV field separator (default: ";")')
@click.option('--decimal-sep', default=',', 
             help='Decimal separator for numeric values (default: ",")')
@click.pass_context
def main(ctx, output_dir, csv_sep, decimal_sep):
    """MHC-I Binding Prediction and Analysis Tool
    
    A comprehensive CLI for peptide-MHC binding prediction, immunogenicity analysis
    and variant generation. Supports multiple prediction methods and allele types.
    """
    ctx.obj = {
        'predictor': IEDBBindingPredictor(
            output_dir=output_dir,
            csv_separator=csv_sep,
            decimal_separator=decimal_sep
        )
    }

@main.command(epilog="""
EXAMPLES:
  From file:
  mhci_predict predict --peptides peptides.txt --alleles HLA-A*02:01 --method netmhcpan_ba
  
  Direct input:
  mhci_predict predict --peptides "SIINFEKL,RAKFKQLL" --alleles HLA-B*07:02
""")
@click.option('--peptides', 
             help='Input peptides (file path or comma-separated list)')
@click.option('--alleles', required=True,
             help='Comma-separated list of MHC alleles (e.g. "HLA-A*02:01,HLA-B*07:02")')
@click.option('--length', default=9, show_default=True,
             help='Peptide length (8-15 residues)')
@click.option('--method', default='netmhcpan_el', show_default=True,
             help="Prediction method: netmhcpan_el (eluted ligand) or netmhcpan_ba (binding affinity)")
@click.pass_context
def predict(ctx, peptides, alleles, length, method):
    """Run MHC binding predictions
    
    Predict peptide-MHC binding using either eluted ligand (EL) or
    binding affinity (BA) methods. Supports batch processing of peptides.
    """
    predictor = ctx.obj['predictor']
    
    # Read peptides
    if peptides and os.path.exists(peptides):
        with open(peptides) as f:
            peptides_list = [line.strip() for line in f]
    else:
        peptides_list = [p.strip() for p in peptides.split(',')]
    
    # Process alleles
    alleles_list = alleles.split(',')
    
    results = predictor.analyze_peptides_batch(
        peptides_list,
        alleles_list,
        length=length
    )
    
    click.echo(f"Results saved to {predictor.output_dir}")

@main.command(epilog="""
EXAMPLES:
  Filter by multiple criteria:
  mhci_predict filter --input results.csv --score 0.7 --ic50 500 --immunogenicity 1.0
  
  Save filtered results:
  mhci_predict filter --input results.csv --rank 2.0 --output filtered.csv
""")
@click.option('--input', required=True,
             help='Input CSV file from predictions')
@click.option('--output',
             help='Output file path (default: print to console)')
@click.option('--score', type=float,
             help='Minimum binding score (higher = better binding)')
@click.option('--rank', type=float,
             help='Maximum percentile rank (lower = better binding)')
@click.option('--ic50', type=float,
             help='Maximum IC50 nM value (lower = better binding)')
@click.option('--immunogenicity', type=float,
             help='Minimum immunogenicity score (higher = more immunogenic)')
@click.pass_context
def filter(ctx, input, output, score, rank, ic50, immunogenicity):
    """Filter prediction results by binding criteria
    
    Apply multiple thresholds to identify strong binders. All specified
    criteria must be satisfied (logical AND). Supports saving to file
    or console output.
    """
    predictor = ctx.obj['predictor']
    df = pd.read_csv(input, sep=predictor.csv_separator)
    
    filtered = predictor.filter_binders(
        df,
        threshold=score,
        rank_threshold=rank,
        ic50_threshold=ic50,
        immunogenicity_threshold=immunogenicity
    )
    
    if output:
        predictor.save_to_csv(filtered, output)
    else:
        click.echo(filtered.to_string())

@main.command(epilog="""
EXAMPLES:
  Basic pattern:
  mhci_predict generate-variants "A[CD]E[FY]GH"
  
  Save to file:
  mhci_predict generate-variants "G[ILV][AG]STV" --output variants.txt
""")
@click.argument('pattern')
@click.option('--output', 
             help='Output file (default: print to console)')
def generate_variants(pattern, output):
    """Generate peptide variants from sequence pattern
    
    PATTERN: Sequence with variant positions in brackets (e.g. A[CD]E[FY]GH)
    - Generates all combinatorial variants
    - Supports single-letter amino acid codes
    """
    predictor = IEDBBindingPredictor()  # Temporary instance
    variants = predictor.generate_variants(pattern)
    
    if output:
        with open(output, 'w') as f:
            f.write('\n'.join(variants))
    else:
        click.echo('\n'.join(variants))

@main.command(epilog="""
EXAMPLES:
  European format:
  mhci_predict set-config --csv-sep ";" --decimal-sep ","
  
  Reset to default:
  mhci_predict set-config --csv-sep "," --decimal-sep "."
""")
@click.option('--csv-sep', 
             help='Set CSV field separator (",", ";")')
@click.option('--decimal-sep',
             help='Set decimal separator (".", ",")')
@click.pass_context
def set_config(ctx, csv_sep, decimal_sep):
    """Configure CSV format settings
    
    Adjust output formatting for regional preferences. Changes apply
    to subsequent commands in the same session.
    """
    predictor = ctx.obj['predictor']
    if csv_sep:
        predictor.set_csv_format(csv_separator=csv_sep)
    if decimal_sep:
        predictor.set_csv_format(decimal_separator=decimal_sep)
    click.echo("Configuration updated")

@main.command(epilog="""
EXAMPLES:
  Using pattern:
  mhci_predict full-analysis --pattern "A[CD]E[FY]GH" --alleles HLA-A*02:01 --output results.csv
  
  Using peptide list:
  mhci_predict full-analysis --peptides peptides.txt --alleles HLA-B*07:02,HLA-A*24:02 --length 10 --output data.csv
""")
@click.option('--peptides', 
             help='Peptides input (file or comma-separated list)')
@click.option('--pattern',
             help='Pattern for generating variants (e.g. "A[CD]E[FY]GH")')
@click.option('--alleles', required=True,
             help='Comma-separated list of MHC alleles')
@click.option('--length', default=9, show_default=True,
             help='Peptide length (required for predictions)')
@click.option('--output', required=True,
             help='Output CSV file path')
@click.pass_context
def full_analysis(ctx, peptides, pattern, alleles, length, output):
    """Complete analysis pipeline
    
    Combines variant generation, binding predictions (EL + BA), and
    immunogenicity scoring in one workflow. Produces comprehensive
    results with all critical metrics.
    """
    predictor = ctx.obj['predictor']
    
    # Generate or read peptides
    if pattern:
        peptides_list = predictor.generate_variants(pattern)
    elif peptides:
        if os.path.exists(peptides):
            with open(peptides) as f:
                peptides_list = [line.strip() for line in f]
        else:
            peptides_list = [p.strip() for p in peptides.split(',')]
    else:
        raise click.UsageError("Must provide either --peptides or --pattern")
    
    alleles_list = alleles.split(',')
    
    # Get predictions for all alleles
    results = predictor.analyze_peptides_batch(peptides_list, alleles_list, length)
    
    # Save the results
    predictor.save_to_csv(results, output)
    click.echo(f"✅ Full analysis completed. Results saved to {output}")

if __name__ == '__main__':
    main()
