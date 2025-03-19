"""
Peptide Analysis Library

Una libreria per la generazione, analisi e visualizzazione di peptidi 9-mer.
"""

__version__ = "0.1.0"

from .peptide_generator import download_swissprot_data, generate_9mers, save_to_csv
from .peptide_analyzer import analyze_peptides, parse_iedb_response
from .visualization import plot_score_distribution, plot_percentile_distribution, plot_allele_comparison, plot_immunogenicity_correlation, plot_category_distribution
from .utils import load_peptides_from_csv, load_results_from_csv, filter_results_by_percentile, filter_results_by_immunogenicity, filter_results_combined
from .sequence_variants import generate_sequence_variants, save_variants_to_csv, generate_all_variants
from .immunogenicity import (
    predict_peptide_immunogenicity, 
    get_available_immunogenicity_alleles,
    print_available_immunogenicity_alleles
)

# Funzione principale che combina tutte le funzionalità
from .main import (
    generate_and_analyze_peptides, 
    run_complete_analysis, 
    analyze_sequence_variants,
    analyze_peptide_immunogenicity
)
