import os
import click
import pandas as pd
from .core import IEDBBindingPredictor

@click.group()
@click.option('--output-dir', default='./output', help='Output directory')
@click.option('--csv-sep', default=',', help='CSV separator')
@click.option('--decimal-sep', default='.', help='Decimal separator')
@click.pass_context
def main(ctx, output_dir, csv_sep, decimal_sep):
    """MHC-I Binding Prediction CLI Tool"""
    ctx.obj = {
        'predictor': IEDBBindingPredictor(
            output_dir=output_dir,
            csv_separator=csv_sep,
            decimal_separator=decimal_sep
        )
    }

@main.command()
@click.option('--peptides', help='Peptide file or comma-separated list')
@click.option('--alleles', required=True, help='Comma-separated list of alleles')
@click.option('--length', default=9, help='Peptide length')
@click.option('--method', default='netmhcpan_el', help='Prediction method')
@click.pass_context
def predict(ctx, peptides, alleles, length, method):
    """Run binding predictions"""
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

@main.command()
@click.option('--input', required=True, help='Input CSV file')
@click.option('--output', help='Output CSV file')
@click.option('--score', type=float, help='Score threshold')
@click.option('--rank', type=float, help='Percentile rank threshold')
@click.option('--ic50', type=float, help='IC50 threshold')
@click.option('--immunogenicity', type=float, help='Immunogenicity threshold')
@click.pass_context
def filter(ctx, input, output, score, rank, ic50, immunogenicity):
    """Filter binding results"""
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

@main.command()
@click.argument('pattern')
@click.option('--output', help='Output file')
def generate_variants(pattern, output):
    """Generate peptide variants from pattern"""
    predictor = IEDBBindingPredictor()  # Temporary instance
    variants = predictor.generate_variants(pattern)
    
    if output:
        with open(output, 'w') as f:
            f.write('\n'.join(variants))
    else:
        click.echo('\n'.join(variants))

@main.command()
@click.option('--csv-sep', help='Set CSV separator')
@click.option('--decimal-sep', help='Set decimal separator')
@click.pass_context
def set_config(ctx, csv_sep, decimal_sep):
    """Update configuration settings"""
    predictor = ctx.obj['predictor']
    if csv_sep:
        predictor.set_csv_format(csv_separator=csv_sep)
    if decimal_sep:
        predictor.set_csv_format(decimal_separator=decimal_sep)
    click.echo("Configuration updated")

@main.command()
@click.option('--peptides', help='Peptide file or comma-separated list')
@click.option('--pattern', help='Pattern for generating variants (e.g. A[CD]E[FY]GH)')
@click.option('--alleles', required=True, help='Comma-separated list of alleles')
@click.option('--length', default=9, help='Peptide length')
@click.option('--output', required=True, help='Output CSV file')
@click.pass_context
def full_analysis(ctx, peptides, pattern, alleles, length, output):
    """Complete analysis pipeline with variants generation, predictions and immunogenicity"""
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
