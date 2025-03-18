#!/usr/bin/env python3
"""
Esempio di utilizzo della libreria peptide_analysis.
"""

import os
import argparse
from peptide_analysis import (
    generate_and_analyze_peptides,
    run_complete_analysis,
    load_results_from_csv,
    filter_results_by_percentile,
    plot_score_distribution,
    plot_percentile_distribution,
    plot_allele_comparison
)

def main():
    parser = argparse.ArgumentParser(description="Analisi di peptidi 9-mer")
    
    # Argomenti generali
    parser.add_argument("--mode", choices=["generate", "analyze", "visualize"], default="analyze",
                        help="Modalità di esecuzione")
    parser.add_argument("--output-dir", default="output",
                        help="Directory di output")
    
    # Argomenti per la generazione
    parser.add_argument("--num-peptides", type=int, default=1000,
                        help="Numero di peptidi da generare")
    
    # Argomenti per l'analisi
    parser.add_argument("--input-csv", 
                        help="File CSV di input con i peptidi")
    parser.add_argument("--alleles", default="HLA-A*01:01,HLA-A*02:01",
                        help="Lista di alleli HLA separati da virgola")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Dimensione del batch per le richieste")
    parser.add_argument("--percentile-threshold", type=float, default=10.0,
                        help="Soglia di percentile rank per il filtraggio")
    parser.add_argument("--percentile-operator", default="<=", choices=["<", "<=", ">", ">=", "=="],
                        help="Operatore di confronto per il filtraggio del percentile rank")
    
    # Argomenti per la visualizzazione
    parser.add_argument("--results-csv", 
                        help="File CSV con i risultati dell'analisi")
    
    args = parser.parse_args()
    
    # Eseguiamo l'operazione richiesta
    if args.mode == "generate":
        # Generiamo e analizziamo i peptidi
        peptides, results = generate_and_analyze_peptides(
            args.num_peptides, args.alleles, args.batch_size, args.output_dir
        )
        print(f"Generati e analizzati {len(peptides)} peptidi")
        
    elif args.mode == "analyze":
        # Eseguiamo l'analisi completa
        report = run_complete_analysis(
            args.input_csv, args.num_peptides, args.alleles, 
            args.batch_size, args.output_dir, args.percentile_threshold,
            args.percentile_operator
        )
        print("\nReport dell'analisi:")
        for key, value in report.items():
            print(f"- {key}: {value}")
        
    elif args.mode == "visualize":
        # Visualizziamo i risultati
        if not args.results_csv:
            print("Errore: è necessario specificare --results-csv per la modalità visualize")
            return
        
        results = load_results_from_csv(args.results_csv)
        if not results:
            print(f"Nessun risultato trovato in {args.results_csv}")
            return
            
        # Filtriamo i risultati se è stata specificata una soglia
        if args.percentile_threshold is not None:
            results = filter_results_by_percentile(
                results, args.percentile_threshold, args.percentile_operator
            )
        
        plots_dir = os.path.join(args.output_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        plot_score_distribution(results, os.path.join(plots_dir, "score_distribution.png"))
        plot_percentile_distribution(results, os.path.join(plots_dir, "percentile_distribution.png"))
        plot_allele_comparison(results, os.path.join(plots_dir, "allele_comparison.png"))
        
        print(f"Grafici salvati in {plots_dir}")

if __name__ == "__main__":
    main()
