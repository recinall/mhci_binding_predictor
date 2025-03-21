#!/usr/bin/env python
"""
Command Line Interface (CLI) per la libreria peptide_analysis.

Questo modulo fornisce un'interfaccia a riga di comando per le principali
funzionalità della libreria peptide_analysis.
"""

import os
import sys
import argparse
import pandas as pd
import csv
import logging
from datetime import datetime

from .sequence_variants import generate_all_variants
from .peptide_analyzer import analyze_peptides
from .binding_prediction import predict_binding
from .immunogenicity import predict_peptide_immunogenicity
from .utils import (
    filter_results_by_percentile,
    filter_results_by_immunogenicity,
    filter_results_by_ic50,
    filter_results_combined,
    ensure_directory
)
from .visualization import (
    plot_score_distribution,
    plot_percentile_distribution,
    plot_allele_comparison,
    plot_immunogenicity_correlation,
    plot_category_distribution
)

# Configurazione del logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def setup_parser():
    """Configura il parser degli argomenti da linea di comando."""
    parser = argparse.ArgumentParser(
        description='Command Line Interface per peptide_analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Esempi di utilizzo:
  # Generare varianti di sequenze da un file
  peptide-analysis variants --input sequences.txt --output variants.csv
  
  # Analizzare peptidi con IEDB
  peptide-analysis analyze --input peptides.csv --alleles HLA-A*02:01 --output results.csv
  
  # Predire l'immunogenicità
  peptide-analysis immunogenicity --input peptides.csv --allele HLA-A0201 --output immuno_results.csv
  
  # Predire il binding locale con predict_binding.py
  peptide-analysis binding --input peptides.csv --allele HLA-A*02:01 --method netmhcpan_el --output binding_results.csv
  
  # Filtrare risultati
  peptide-analysis filter --input results.csv --percentile 0.5 --immunogenicity 0 --ic50 500 --output filtered.csv
  
  # Generare grafici
  peptide-analysis visualize --input results.csv --output-dir plots
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Comando da eseguire')
    
    # Parser per il comando 'variants'
    variants_parser = subparsers.add_parser('variants', help='Genera varianti di sequenze')
    variants_parser.add_argument('--input', required=True, help='File di input con le sequenze con varianti')
    variants_parser.add_argument('--output', required=True, help='File CSV di output per le varianti generate')
    variants_parser.add_argument('--output-dir', default='variants', help='Directory di output per i file CSV delle singole sequenze')
    variants_parser.add_argument('--csv-separator', default=',', help='Separatore per il file CSV (default: ,)')
    variants_parser.add_argument('--csv-decimal', default='.', help='Carattere decimale per il file CSV (default: .)')
    
    # Parser per il comando 'analyze'
    analyze_parser = subparsers.add_parser('analyze', help='Analizza peptidi con IEDB')
    analyze_parser.add_argument('--input', required=True, help='File CSV di input con i peptidi')
    analyze_parser.add_argument('--alleles', required=True, help='Lista di alleli separati da virgola (es. HLA-A*02:01,HLA-B*07:02)')
    analyze_parser.add_argument('--output', required=True, help='File CSV di output per i risultati')
    analyze_parser.add_argument('--batch-size', type=int, default=10, help='Dimensione del batch per le richieste IEDB (default: 10)')
    analyze_parser.add_argument('--csv-separator', default=',', help='Separatore per il file CSV (default: ,)')
    analyze_parser.add_argument('--csv-decimal', default='.', help='Carattere decimale per il file CSV (default: .)')
    
    # Parser per il comando 'immunogenicity'
    immuno_parser = subparsers.add_parser('immunogenicity', help='Predici immunogenicità dei peptidi')
    immuno_parser.add_argument('--input', required=True, help='File CSV di input con i peptidi')
    immuno_parser.add_argument('--allele', help='Allele per la predizione (es. HLA-A0201)')
    immuno_parser.add_argument('--output', required=True, help='File CSV di output per i risultati')
    immuno_parser.add_argument('--csv-separator', default=',', help='Separatore per il file CSV (default: ,)')
    immuno_parser.add_argument('--csv-decimal', default='.', help='Carattere decimale per il file CSV (default: .)')
    immuno_parser.add_argument('--list-alleles', action='store_true', help='Mostra gli alleli disponibili per la predizione dell\'immunogenicità')
    
    # Parser per il comando 'binding'
    binding_parser = subparsers.add_parser('binding', help='Predici binding locale con predict_binding.py')
    binding_parser.add_argument('--input', required=True, help='File CSV di input con i peptidi')
    binding_parser.add_argument('--allele', required=True, help='Allele per la predizione (es. HLA-A*02:01)')
    binding_parser.add_argument('--method', default='netmhcpan_el', help='Metodo di predizione (default: netmhcpan_el)')
    binding_parser.add_argument('--output', required=True, help='File CSV di output per i risultati')
    binding_parser.add_argument('--predict-binding-path', help='Percorso al file predict_binding.py')
    binding_parser.add_argument('--csv-separator', default=',', help='Separatore per il file CSV (default: ,)')
    binding_parser.add_argument('--csv-decimal', default='.', help='Carattere decimale per il file CSV (default: .)')
    binding_parser.add_argument('--list-methods', action='store_true', help='Mostra i metodi disponibili per la predizione del binding')
    binding_parser.add_argument('--list-alleles', action='store_true', help='Mostra gli alleli disponibili per il metodo specificato')
    
    # Parser per il comando 'filter'
    filter_parser = subparsers.add_parser('filter', help='Filtra risultati')
    filter_parser.add_argument('--input', required=True, help='File CSV di input con i risultati')
    filter_parser.add_argument('--output', required=True, help='File CSV di output per i risultati filtrati')
    filter_parser.add_argument('--percentile', type=float, help='Soglia di percentile rank')
    filter_parser.add_argument('--percentile-operator', default='<', choices=['<', '<=', '>', '>=', '=='], help='Operatore per il percentile rank (default: <)')
    filter_parser.add_argument('--immunogenicity', type=float, help='Soglia di immunogenicità')
    filter_parser.add_argument('--immunogenicity-operator', default='>', choices=['<', '<=', '>', '>=', '=='], help='Operatore per l\'immunogenicità (default: >)')
    filter_parser.add_argument('--ic50', type=float, help='Soglia di IC50')
    filter_parser.add_argument('--ic50-operator', default='<', choices=['<', '<=', '>', '>=', '=='], help='Operatore per l\'IC50 (default: <)')
    filter_parser.add_argument('--csv-separator', default=',', help='Separatore per il file CSV (default: ,)')
    filter_parser.add_argument('--csv-decimal', default='.', help='Carattere decimale per il file CSV (default: .)')
    
    # Parser per il comando 'visualize'
    visualize_parser = subparsers.add_parser('visualize', help='Genera grafici di analisi')
    visualize_parser.add_argument('--input', required=True, help='File CSV di input con i risultati')
    visualize_parser.add_argument('--output-dir', required=True, help='Directory di output per i grafici')
    visualize_parser.add_argument('--csv-separator', default=',', help='Separatore per il file CSV di input (default: ,)')
    visualize_parser.add_argument('--csv-decimal', default='.', help='Carattere decimale per il file CSV di input (default: .)')
    visualize_parser.add_argument('--format', default='png', choices=['png', 'pdf', 'svg', 'jpg'], help='Formato dei grafici (default: png)')
    visualize_parser.add_argument('--dpi', type=int, default=300, help='DPI per i grafici (default: 300)')
    
    # Parser per il comando 'complete'
    complete_parser = subparsers.add_parser('complete', help='Esegui analisi completa (varianti, binding, immunogenicità, filtri, grafici)')
    complete_parser.add_argument('--input', required=True, help='File di input con le sequenze con varianti o lista di peptidi')
    complete_parser.add_argument('--is-variants', action='store_true', help='Indica se il file di input contiene sequenze con varianti')
    complete_parser.add_argument('--allele', required=True, help='Allele per la predizione (es. HLA-A*02:01)')
    complete_parser.add_argument('--method', default='netmhcpan_el', help='Metodo di predizione del binding (default: netmhcpan_el)')
    complete_parser.add_argument('--predict-binding-path', help='Percorso al file predict_binding.py')
    complete_parser.add_argument('--output-dir', required=True, help='Directory di output per i risultati e i grafici')
    complete_parser.add_argument('--percentile', type=float, default=0.5, help='Soglia di percentile rank per il filtraggio (default: 0.5)')
    complete_parser.add_argument('--immunogenicity', type=float, default=0.0, help='Soglia di immunogenicità per il filtraggio (default: 0.0)')
    complete_parser.add_argument('--ic50', type=float, default=500, help='Soglia di IC50 per il filtraggio (default: 500)')
    complete_parser.add_argument('--csv-separator', default=',', help='Separatore per i file CSV (default: ,)')
    complete_parser.add_argument('--csv-decimal', default='.', help='Carattere decimale per i file CSV (default: .)')
    
    return parser

def parse_variants_file(file_path):
    """
    Analizza un file di sequenze con varianti.
    
    Args:
        file_path: Percorso al file di input
        
    Returns:
        list: Lista di sequenze con varianti
    """
    sequences = []
    current_sequence = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = []
            else:
                # Formato atteso: [A,B,C] o [A, B, C] o A,B,C o A, B, C
                line = line.strip('[]')
                variants = [v.strip() for v in line.split(',')]
                current_sequence.append(variants)
    
    # Aggiungi l'ultima sequenza se presente
    if current_sequence:
        sequences.append(current_sequence)
    
    return sequences

def load_peptides_from_csv(file_path, separator=','):
    """
    Carica i peptidi da un file CSV.
    
    Args:
        file_path: Percorso al file CSV
        separator: Separatore del CSV
        
    Returns:
        list: Lista di peptidi
    """
    peptides = []
    
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter=separator)
        header = next(reader, None)
        
        # Cerca la colonna 'peptide' o usa la prima colonna
        peptide_col = 0
        if header:
            for i, col in enumerate(header):
                if col.lower() == 'peptide':
                    peptide_col = i
                    break
        
        for row in reader:
            if row and len(row) > peptide_col:
                peptide = row[peptide_col].strip()
                if peptide:
                    peptides.append(peptide)
    
    return peptides

def save_to_csv(data, file_path, separator=',', decimal='.'):
    """
    Salva i dati in un file CSV.
    
    Args:
        data: Dati da salvare (lista di dizionari o DataFrame)
        file_path: Percorso al file CSV
        separator: Separatore del CSV
        decimal: Carattere decimale
    """
    if isinstance(data, list) and data and isinstance(data[0], dict):
        # Converti la lista di dizionari in DataFrame
        df = pd.DataFrame(data)
    elif isinstance(data, pd.DataFrame):
        df = data
    else:
        logger.error("Formato dati non supportato per il salvataggio in CSV")
        return
    
    # Sostituisci il separatore decimale se necessario
    if decimal != '.':
        for col in df.select_dtypes(include=['float']).columns:
            df[col] = df[col].astype(str).str.replace('.', decimal)
    
    # Salva il DataFrame in CSV
    df.to_csv(file_path, sep=separator, index=False)
    logger.info(f"Dati salvati in {file_path}")

def handle_variants_command(args):
    """Gestisce il comando 'variants'."""
    logger.info(f"Generazione di varianti di sequenze da {args.input}")
    
    # Analizza il file di sequenze con varianti
    sequences = parse_variants_file(args.input)
    
    if not sequences:
        logger.error("Nessuna sequenza trovata nel file di input")
        return
    
    logger.info(f"Trovate {len(sequences)} sequenze con varianti")
    
    # Genera tutte le varianti
    ensure_directory(args.output_dir)
    all_variants = generate_all_variants(sequences, args.output_dir, args.output)
    
    logger.info(f"Generate {len(all_variants)} varianti totali")
    logger.info(f"Varianti salvate in {args.output}")

def handle_analyze_command(args):
    """Gestisce il comando 'analyze'."""
    logger.info(f"Analisi di peptidi da {args.input} con alleli {args.alleles}")
    
    # Carica i peptidi dal file CSV
    peptides = load_peptides_from_csv(args.input, args.csv_separator)
    
    if not peptides:
        logger.error("Nessun peptide trovato nel file di input")
        return
    
    logger.info(f"Trovati {len(peptides)} peptidi")
    
    # Analizza i peptidi con IEDB
    results = analyze_peptides(
        peptides=peptides,
        allele_list=args.alleles,
        batch_size=args.batch_size,
        output_csv=args.output
    )
    
    logger.info(f"Analisi completata. Risultati salvati in {args.output}")
    
    return results

def handle_immunogenicity_command(args):
    """Gestisce il comando 'immunogenicity'."""
    if args.list_alleles:
        from .immunogenicity import print_available_immunogenicity_alleles
        print_available_immunogenicity_alleles()
        return
    
    logger.info(f"Predizione dell'immunogenicità per i peptidi in {args.input}")
    
    # Carica i peptidi dal file CSV
    peptides = load_peptides_from_csv(args.input, args.csv_separator)
    
    if not peptides:
        logger.error("Nessun peptide trovato nel file di input")
        return
    
    logger.info(f"Trovati {len(peptides)} peptidi")
    
    # Predici l'immunogenicità
    results = predict_peptide_immunogenicity(
        peptides=peptides,
        allele=args.allele,
        output_csv=args.output
    )
    
    logger.info(f"Predizione completata. Risultati salvati in {args.output}")
    
    return results

def handle_binding_command(args):
    """Gestisce il comando 'binding'."""
    if args.list_methods:
        from .binding_prediction import get_available_methods
        methods = get_available_methods(args.predict_binding_path)
        print("Metodi disponibili per la predizione del binding:")
        for method in methods:
            print(f"- {method}")
        return
    
    if args.list_alleles:
        from .binding_prediction import get_available_alleles
        alleles = get_available_alleles(args.method, args.predict_binding_path)
        print(f"Alleli disponibili per il metodo {args.method}:")
        print(alleles)
        return
    
    logger.info(f"Predizione del binding per i peptidi in {args.input} con allele {args.allele} e metodo {args.method}")
    
    # Carica i peptidi dal file CSV
    peptides = load_peptides_from_csv(args.input, args.csv_separator)
    
    if not peptides:
        logger.error("Nessun peptide trovato nel file di input")
        return
    
    logger.info(f"Trovati {len(peptides)} peptidi")
    
    # Raggruppa i peptidi per lunghezza
    peptides_by_length = {}
    for peptide in peptides:
        length = len(peptide)
        if length not in peptides_by_length:
            peptides_by_length[length] = []
        peptides_by_length[length].append(peptide)
    
    # Predici il binding per ciascun gruppo di peptidi
    all_results = pd.DataFrame()
    
    for length, peptides_group in peptides_by_length.items():
        logger.info(f"Predizione per {len(peptides_group)} peptidi di lunghezza {length}")
        
        try:
            results = predict_binding(
                peptides=peptides_group,
                allele=args.allele,
                method=args.method,
                predict_binding_path=args.predict_binding_path
            )
            
            if not results.empty:
                all_results = pd.concat([all_results, results])
        except Exception as e:
            logger.error(f"Errore nella predizione per peptidi di lunghezza {length}: {e}")
    
    # Salva i risultati
    if not all_results.empty:
        save_to_csv(all_results, args.output, args.csv_separator, args.csv_decimal)
        logger.info(f"Predizione completata. Risultati salvati in {args.output}")
    else:
        logger.error("Nessun risultato ottenuto dalla predizione")
    
    return all_results

def handle_filter_command(args):
    """Gestisce il comando 'filter'."""
    logger.info(f"Filtraggio dei risultati da {args.input}")
    
    # Carica i risultati dal file CSV
    try:
        df = pd.read_csv(args.input, sep=args.csv_separator)
        
        # Converti il separatore decimale se necessario
        if args.csv_decimal != '.':
            for col in df.select_dtypes(include=['object']).columns:
                try:
                    df[col] = df[col].str.replace(args.csv_decimal, '.').astype(float)
                except:
                    pass
        
        # Converti il DataFrame in lista di dizionari
        results = df.to_dict('records')
        
        logger.info(f"Caricati {len(results)} risultati")
    except Exception as e:
        logger.error(f"Errore nel caricamento dei risultati: {e}")
        return
    
    # Applica i filtri
    filtered_results = results
    
    if args.percentile is not None:
        filtered_results = filter_results_by_percentile(
            filtered_results,
            threshold=args.percentile,
            operator=args.percentile_operator
        )
    
    if args.immunogenicity is not None:
        filtered_results = filter_results_by_immunogenicity(
            filtered_results,
            threshold=args.immunogenicity,
            operator=args.immunogenicity_operator
        )
    
    if args.ic50 is not None:
        filtered_results = filter_results_by_ic50(
            filtered_results,
            threshold=args.ic50,
            operator=args.ic50_operator
        )
    
    # Salva i risultati filtrati
    save_to_csv(filtered_results, args.output, args.csv_separator, args.csv_decimal)
    logger.info(f"Filtraggio completato. Risultati salvati in {args.output}")
    
    return filtered_results

def handle_visualize_command(args):
    """Gestisce il comando 'visualize'."""
    logger.info(f"Generazione di grafici per i risultati in {args.input}")
    
    # Carica i risultati dal file CSV
    try:
        df = pd.read_csv(args.input, sep=args.csv_separator)
        
        # Converti il separatore decimale se necessario
        if args.csv_decimal != '.':
            for col in df.select_dtypes(include=['object']).columns:
                try:
                    df[col] = df[col].str.replace(args.csv_decimal, '.').astype(float)
                except:
                    pass
        
        # Converti il DataFrame in lista di dizionari
        results = df.to_dict('records')
        
        logger.info(f"Caricati {len(results)} risultati")
    except Exception as e:
        logger.error(f"Errore nel caricamento dei risultati: {e}")
        return
    
    # Crea la directory di output
    ensure_directory(args.output_dir)
    
    # Genera i grafici
    file_extension = f".{args.format}"
    
    try:
        # Grafico della distribuzione degli score
        score_plot = os.path.join(args.output_dir, f"score_distribution{file_extension}")
        plot_score_distribution(results, score_plot, dpi=args.dpi)
        logger.info(f"Generato grafico della distribuzione degli score: {score_plot}")
        
        # Grafico della distribuzione dei percentile rank
        percentile_plot = os.path.join(args.output_dir, f"percentile_distribution{file_extension}")
        plot_percentile_distribution(results, percentile_plot, dpi=args.dpi)
        logger.info(f"Generato grafico della distribuzione dei percentile rank: {percentile_plot}")
        
        # Grafico di confronto tra alleli
        allele_plot = os.path.join(args.output_dir, f"allele_comparison{file_extension}")
        plot_allele_comparison(results, allele_plot, dpi=args.dpi)
        logger.info(f"Generato grafico di confronto tra alleli: {allele_plot}")
        
        # Grafico della correlazione tra immunogenicità e binding
        immuno_plot = os.path.join(args.output_dir, f"immunogenicity_correlation{file_extension}")
        plot_immunogenicity_correlation(results, immuno_plot, dpi=args.dpi)
        logger.info(f"Generato grafico della correlazione tra immunogenicità e binding: {immuno_plot}")
        
        # Grafico della distribuzione delle categorie
        category_plot = os.path.join(args.output_dir, f"category_distribution{file_extension}")
        plot_category_distribution(results, category_plot, dpi=args.dpi)
        logger.info(f"Generato grafico della distribuzione delle categorie: {category_plot}")
    
    except Exception as e:
        logger.error(f"Errore nella generazione dei grafici: {e}")

def handle_complete_command(args):
    """Gestisce il comando 'complete'."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(args.output_dir, f"analysis_{timestamp}")
    ensure_directory(output_dir)
    
    logger.info(f"Esecuzione dell'analisi completa con output in {output_dir}")
    
    # Step 1: Genera varianti o carica peptidi
    if args.is_variants:
        logger.info("Generazione di varianti di sequenze")
        sequences = parse_variants_file(args.input)
        
        if not sequences:
            logger.error("Nessuna sequenza trovata nel file di input")
            return
        
        variants_file = os.path.join(output_dir, "variants.csv")
        variants_dir = os.path.join(output_dir, "variants")
        ensure_directory(variants_dir)
        
        all_variants = generate_all_variants(sequences, variants_dir, variants_file)
        peptides = all_variants
        
        logger.info(f"Generate {len(peptides)} varianti totali")
    else:
        logger.info("Caricamento di peptidi dal file CSV")
        peptides = load_peptides_from_csv(args.input, args.csv_separator)
        
        if not peptides:
            logger.error("Nessun peptide trovato nel file di input")
            return
        
        logger.info(f"Caricati {len(peptides)} peptidi")
    
    # Step 2: Predici il binding
    binding_file = os.path.join(output_dir, "binding_results.csv")
    
    # Raggruppa i peptidi per lunghezza
    peptides_by_length = {}
    for peptide in peptides:
        length = len(peptide)
        if length not in peptides_by_length:
            peptides_by_length[length] = []
        peptides_by_length[length].append(peptide)
    
    # Predici il binding per ciascun gruppo di peptidi
    all_binding_results = pd.DataFrame()
    
    for length, peptides_group in peptides_by_length.items():
        logger.info(f"Predizione del binding per {len(peptides_group)} peptidi di lunghezza {length}")
        
        try:
            binding_results = predict_binding(
                peptides=peptides_group,
                allele=args.allele,
                method=args.method,
                predict_binding_path=args.predict_binding_path
            )
            
            if not binding_results.empty:
                all_binding_results = pd.concat([all_binding_results, binding_results])
        except Exception as e:
            logger.error(f"Errore nella predizione del binding per peptidi di lunghezza {length}: {e}")
    
    # Salva i risultati del binding
    if not all_binding_results.empty:
        save_to_csv(all_binding_results, binding_file, args.csv_separator, args.csv_decimal)
        logger.info(f"Risultati del binding salvati in {binding_file}")
    else:
        logger.error("Nessun risultato ottenuto dalla predizione del binding")
        return
    
    # Step 3: Predici l'immunogenicità
    immuno_file = os.path.join(output_dir, "immunogenicity_results.csv")
    
    logger.info("Predizione dell'immunogenicità")
    
    # Estrai i peptidi dai risultati del binding
    binding_peptides = all_binding_results['peptide'].tolist() if 'peptide' in all_binding_results.columns else []
    
    if binding_peptides:
        try:
            immuno_results = predict_peptide_immunogenicity(
                peptides=binding_peptides,
                allele=args.allele.replace('*', ''),  # Adatta il formato dell'allele per l'immunogenicità
                output_csv=immuno_file
            )
            
            logger.info(f"Risultati dell'immunogenicità salvati in {immuno_file}")
        except Exception as e:
            logger.error(f"Errore nella predizione dell'immunogenicità: {e}")
            immuno_results = []
    else:
        logger.error("Nessun peptide disponibile per la predizione dell'immunogenicità")
        immuno_results = []
    
    # Step 4: Combina i risultati
    combined_file = os.path.join(output_dir, "combined_results.csv")
    
    logger.info("Combinazione dei risultati di binding e immunogenicità")
    
    # Crea un dizionario per i risultati di immunogenicità
    immuno_dict = {}
    for result in immuno_results:
        immuno_dict[result['peptide']] = result['score']
    
    # Combina i risultati
    combined_results = []
    
    for i, row in all_binding_results.iterrows():
        if 'peptide' not in row:
            continue
        
        peptide = row['peptide']
        
        result = {
            'peptide': peptide,
            'allele': args.allele,
            'percentile_rank': row['rank'] if 'rank' in row else None,
            'score': row['score'] if 'score' in row else None
        }
        
        # Aggiungi IC50 se disponibile
        if 'ic50' in row:
            result['ic50'] = row['ic50']
        elif 'ann_ic50' in row:
            result['ic50'] = row['ann_ic50']
        
        # Aggiungi immunogenicity score
        if peptide in immuno_dict:
            result['immunogenicity_score'] = immuno_dict[peptide]
        
        combined_results.append(result)
    
    # Salva i risultati combinati
    save_to_csv(combined_results, combined_file, args.csv_separator, args.csv_decimal)
    logger.info(f"Risultati combinati salvati in {combined_file}")
    
    # Step 5: Filtra i risultati
    filtered_file = os.path.join(output_dir, "filtered_results.csv")
    
    logger.info("Filtraggio dei risultati")
    
    # Applica i filtri
    filtered_results = combined_results
    
    filtered_results = filter_results_by_percentile(
        filtered_results,
        threshold=args.percentile,
        operator="<"
    )
    
    filtered_results = filter_results_by_immunogenicity(
        filtered_results,
        threshold=args.immunogenicity,
        operator=">"
    )
    
    filtered_results = filter_results_by_ic50(
        filtered_results,
        threshold=args.ic50,
        operator="<"
    )
    
    # Salva i risultati filtrati
    save_to_csv(filtered_results, filtered_file, args.csv_separator, args.csv_decimal)
    logger.info(f"Risultati filtrati salvati in {filtered_file}")
    
    # Step 6: Genera grafici
    plots_dir = os.path.join(output_dir, "plots")
    ensure_directory(plots_dir)
    
    logger.info("Generazione dei grafici")
    
    try:
        # Grafico della distribuzione degli score
        score_plot = os.path.join(plots_dir, "score_distribution.png")
        plot_score_distribution(combined_results, score_plot)
        logger.info(f"Generato grafico della distribuzione degli score: {score_plot}")
        
        # Grafico della distribuzione dei percentile rank
        percentile_plot = os.path.join(plots_dir, "percentile_distribution.png")
        plot_percentile_distribution(combined_results, percentile_plot)
        logger.info(f"Generato grafico della distribuzione dei percentile rank: {percentile_plot}")
        
        # Grafico di confronto tra alleli
        allele_plot = os.path.join(plots_dir, "allele_comparison.png")
        plot_allele_comparison(combined_results, allele_plot)
        logger.info(f"Generato grafico di confronto tra alleli: {allele_plot}")
        
        # Grafico della correlazione tra immunogenicità e binding
        immuno_plot = os.path.join(plots_dir, "immunogenicity_correlation.png")
        plot_immunogenicity_correlation(combined_results, immuno_plot)
        logger.info(f"Generato grafico della correlazione tra immunogenicità e binding: {immuno_plot}")
        
        # Grafico della distribuzione delle categorie
        category_plot = os.path.join(plots_dir, "category_distribution.png")
        plot_category_distribution(combined_results, category_plot)
        logger.info(f"Generato grafico della distribuzione delle categorie: {category_plot}")
    
    except Exception as e:
        logger.error(f"Errore nella generazione dei grafici: {e}")
    
    logger.info(f"Analisi completa terminata. Tutti i risultati sono disponibili in {output_dir}")

def main():
    """Funzione principale per la CLI."""
    parser = setup_parser()
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return
    
    # Esegui il comando appropriato
    if args.command == 'variants':
        handle_variants_command(args)
    elif args.command == 'analyze':
        handle_analyze_command(args)
    elif args.command == 'immunogenicity':
        handle_immunogenicity_command(args)
    elif args.command == 'binding':
        handle_binding_command(args)
    elif args.command == 'filter':
        handle_filter_command(args)
    elif args.command == 'visualize':
        handle_visualize_command(args)
    elif args.command == 'complete':
        handle_complete_command(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
