"""
Modulo principale che combina tutte le funzionalità della libreria.
"""

import os
import csv
from .peptide_generator import generate_peptides_pipeline
from .peptide_analyzer import analyze_peptides, process_peptides_from_csv
from .visualization import plot_score_distribution, plot_percentile_distribution, plot_allele_comparison
from .utils import load_peptides_from_csv, load_results_from_csv, ensure_directory, filter_results_by_percentile
from .sequence_variants import generate_sequence_variants, save_variants_to_csv, generate_all_variants
from .immunogenicity import predict_peptide_immunogenicity, get_available_immunogenicity_alleles

def generate_and_analyze_peptides(num_peptides=1000, allele_list="HLA-A*01:01,HLA-A*02:01", 
                                 batch_size=10, output_dir="output"):
    """
    Pipeline completa per generare e analizzare peptidi.
    
    Parametri:
    num_peptides (int): Numero di peptidi da generare
    allele_list (str): Lista di alleli HLA da utilizzare
    batch_size (int): Dimensione del batch per le richieste
    output_dir (str): Directory di output
    
    Returns:
    tuple: (peptidi generati, risultati dell'analisi)
    """
    # Assicuriamo che la directory di output esista
    ensure_directory(output_dir)
    
    # Generiamo i peptidi
    peptides_csv = os.path.join(output_dir, "peptides.csv")
    peptides, _ = generate_peptides_pipeline(num_peptides, peptides_csv)
    
    # Analizziamo i peptidi
    results_csv = os.path.join(output_dir, "results.csv")
    results = analyze_peptides(peptides, allele_list, batch_size, results_csv)
    
    return peptides, results

def run_complete_analysis(input_csv=None, num_peptides=1000, allele_list="HLA-A*01:01,HLA-A*02:01", 
                         batch_size=10, output_dir="output", percentile_threshold=10.0, 
                         percentile_operator="<=", include_immunogenicity=True):
    """
    Esegue un'analisi completa dei peptidi, inclusa la visualizzazione.
    
    Parametri:
    input_csv (str): File CSV di input (opzionale)
    num_peptides (int): Numero di peptidi da generare (se input_csv non è specificato)
    allele_list (str): Lista di alleli HLA da utilizzare
    batch_size (int): Dimensione del batch per le richieste
    output_dir (str): Directory di output
    percentile_threshold (float): Soglia di percentile rank per il filtraggio
    
    Returns:
    dict: Dizionario con i risultati dell'analisi
    """
    # Assicuriamo che la directory di output esista
    ensure_directory(output_dir)
    
    # Determiniamo se dobbiamo generare peptidi o caricarli da un file
    if input_csv and os.path.exists(input_csv):
        print(f"Caricamento peptidi da {input_csv}")
        peptides = load_peptides_from_csv(input_csv)
        results_csv = os.path.join(output_dir, "results.csv")
        
        # Analizziamo i peptidi
        results = analyze_peptides(peptides, allele_list, batch_size, results_csv)
    else:
        print(f"Generazione di {num_peptides} peptidi")
        peptides, results = generate_and_analyze_peptides(
            num_peptides, allele_list, batch_size, output_dir
        )
    
    # Filtriamo i risultati
    filtered_results = filter_results_by_percentile(results, percentile_threshold, percentile_operator)
    
    # Aggiungiamo l'analisi di immunogenicità se richiesto
    if include_immunogenicity:
        print("\nAggiunta dell'analisi di immunogenicità...")
        # Estrai i peptidi unici dai risultati
        unique_peptides = list(set([r['peptide'] for r in results]))
        
        # Directory per i risultati di immunogenicità
        immuno_dir = os.path.join(output_dir, "immunogenicity")
        ensure_directory(immuno_dir)
        
        # Analizziamo l'immunogenicità dei peptidi
        immuno_results = predict_peptide_immunogenicity(
            unique_peptides, 
            output_csv=os.path.join(immuno_dir, "immunogenicity_results.csv")
        )
        
        # Creiamo un dizionario per mappare peptidi a score di immunogenicità
        immuno_scores = {r['peptide']: r['score'] for r in immuno_results}
        
        # Aggiungiamo lo score di immunogenicità ai risultati originali
        for r in results:
            r['immunogenicity_score'] = immuno_scores.get(r['peptide'], None)
        
        # Aggiungiamo lo score di immunogenicità ai risultati filtrati
        for r in filtered_results:
            r['immunogenicity_score'] = immuno_scores.get(r['peptide'], None)
        
        # Salviamo i risultati completi in un nuovo file CSV
        complete_results_csv = os.path.join(output_dir, "complete_results.csv")
        with open(complete_results_csv, 'w', newline='') as csvfile:
            fieldnames = ['peptide', 'allele', 'score', 'percentile_rank', 'immunogenicity_score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in results:
                writer.writerow(result)
        
        print(f"Risultati completi con immunogenicità salvati in {complete_results_csv}")
    
    # Creiamo i grafici
    plots_dir = os.path.join(output_dir, "plots")
    ensure_directory(plots_dir)
    
    score_plot = plot_score_distribution(
        results, os.path.join(plots_dir, "score_distribution.png")
    )
    
    percentile_plot = plot_percentile_distribution(
        results, os.path.join(plots_dir, "percentile_distribution.png")
    )
    
    allele_plot = plot_allele_comparison(
        results, os.path.join(plots_dir, "allele_comparison.png")
    )
    
    # Aggiungiamo il grafico di correlazione tra binding e immunogenicità se abbiamo i dati
    if include_immunogenicity:
        immuno_plot = plot_immunogenicity_correlation(
            results, os.path.join(plots_dir, "immunogenicity_correlation.png")
        )
    
    # Creiamo un report con i risultati
    report = {
        "total_peptides": len(peptides),
        "total_results": len(results),
        "filtered_results": len(filtered_results),
        "alleles": allele_list.split(","),
        "percentile_threshold": percentile_threshold,
        "percentile_operator": percentile_operator,
        "output_directory": os.path.abspath(output_dir),
        "include_immunogenicity": include_immunogenicity
    }
    
    if include_immunogenicity:
        report["immunogenicity_results"] = len(immuno_results)
    
    print("\nAnalisi completata:")
    print(f"- Peptidi analizzati: {report['total_peptides']}")
    print(f"- Risultati totali: {report['total_results']}")
    print(f"- Risultati con percentile rank {percentile_operator} {percentile_threshold}: {report['filtered_results']}")
    if include_immunogenicity:
        print(f"- Risultati di immunogenicità: {report['immunogenicity_results']}")
        print(f"- Risultati completi con immunogenicità: {complete_results_csv}")
    print(f"- Grafici salvati in: {plots_dir}")
    
    return report

def analyze_peptide_immunogenicity(peptides=None, input_csv=None, custom_mask=None, 
                                  allele=None, output_dir="immunogenicity_output"):
    """
    Analizza l'immunogenicità di peptidi.
    
    Parametri:
    peptides (list): Lista di peptidi da analizzare (opzionale)
    input_csv (str): File CSV di input contenente i peptidi (opzionale)
    custom_mask (str): Maschera personalizzata (opzionale)
    allele (str): Allele HLA da utilizzare (opzionale)
    output_dir (str): Directory di output
    
    Returns:
    dict: Dizionario con i risultati dell'analisi
    """
    # Assicuriamo che la directory di output esista
    ensure_directory(output_dir)
    
    # Determiniamo se dobbiamo caricare peptidi da un file
    if not peptides and input_csv and os.path.exists(input_csv):
        print(f"Caricamento peptidi da {input_csv}")
        peptides = load_peptides_from_csv(input_csv)
    
    if not peptides:
        print("Nessun peptide da analizzare.")
        return {"error": "Nessun peptide da analizzare"}
    
    # File di output per i risultati
    output_csv = os.path.join(output_dir, "immunogenicity_results.csv")
    
    # Analizziamo l'immunogenicità dei peptidi
    print(f"Analisi dell'immunogenicità di {len(peptides)} peptidi...")
    results = predict_peptide_immunogenicity(peptides, custom_mask, allele, output_csv)
    
    # Creiamo un report con i risultati
    report = {
        "total_peptides": len(peptides),
        "total_results": len(results),
        "allele": allele,
        "custom_mask": custom_mask,
        "output_directory": os.path.abspath(output_dir),
        "output_csv": output_csv
    }
    
    print("\nAnalisi dell'immunogenicità completata:")
    print(f"- Peptidi analizzati: {report['total_peptides']}")
    print(f"- Risultati totali: {report['total_results']}")
    if allele:
        print(f"- Allele utilizzato: {allele}")
    if custom_mask:
        print(f"- Maschera personalizzata: {custom_mask}")
    print(f"- Risultati salvati in: {output_csv}")
    
    return report

def analyze_sequence_variants(sequences, allele_list="HLA-A*01:01,HLA-A*02:01", 
                             batch_size=10, output_dir="variants_analysis", 
                             percentile_threshold=10.0, percentile_operator="<="):
    """
    Genera varianti di sequenze peptidiche e le analizza.
    
    Parametri:
    sequences (list): Lista di sequenze, dove ogni sequenza è una lista di liste
    allele_list (str): Lista di alleli HLA da utilizzare
    batch_size (int): Dimensione del batch per le richieste
    output_dir (str): Directory di output
    percentile_threshold (float): Soglia di percentile rank per il filtraggio
    percentile_operator (str): Operatore di confronto per il filtraggio
    
    Returns:
    dict: Dizionario con i risultati dell'analisi
    """
    # Assicuriamo che la directory di output esista
    ensure_directory(output_dir)
    
    # Directory per le varianti
    variants_dir = os.path.join(output_dir, "variants")
    ensure_directory(variants_dir)
    
    # Generiamo le varianti
    print("Generazione delle varianti di sequenze...")
    total_variants, csv_files = generate_all_variants(sequences, variants_dir)
    
    # File combinato con tutte le varianti
    combined_csv = os.path.join(variants_dir, "Sequenze_AIO.csv")
    
    # Analizziamo le varianti
    print(f"Analisi di {total_variants} varianti peptidiche...")
    results_csv = os.path.join(output_dir, "variants_results.csv")
    
    # Carichiamo i peptidi dal file combinato
    peptides = load_peptides_from_csv(combined_csv)
    
    # Analizziamo i peptidi
    results = analyze_peptides(peptides, allele_list, batch_size, results_csv)
    
    # Filtriamo i risultati
    filtered_results = filter_results_by_percentile(results, percentile_threshold, percentile_operator)
    
    # Creiamo i grafici
    plots_dir = os.path.join(output_dir, "plots")
    ensure_directory(plots_dir)
    
    score_plot = plot_score_distribution(
        results, os.path.join(plots_dir, "score_distribution.png")
    )
    
    percentile_plot = plot_percentile_distribution(
        results, os.path.join(plots_dir, "percentile_distribution.png")
    )
    
    allele_plot = plot_allele_comparison(
        results, os.path.join(plots_dir, "allele_comparison.png")
    )
    
    # Creiamo un report con i risultati
    report = {
        "total_variants": total_variants,
        "total_results": len(results),
        "filtered_results": len(filtered_results),
        "alleles": allele_list.split(","),
        "percentile_threshold": percentile_threshold,
        "percentile_operator": percentile_operator,
        "output_directory": os.path.abspath(output_dir),
        "csv_files": csv_files
    }
    
    print("\nAnalisi delle varianti completata:")
    print(f"- Varianti generate: {report['total_variants']}")
    print(f"- Risultati totali: {report['total_results']}")
    print(f"- Risultati con percentile rank {percentile_operator} {percentile_threshold}: {report['filtered_results']}")
    print(f"- Grafici salvati in: {plots_dir}")
    
    return report
