#!/usr/bin/env python3
"""
Esempio di utilizzo della funzionalità di filtraggio dei risultati.
"""

import os
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt

# Aggiungiamo la directory principale al path per importare il pacchetto
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from peptide_analysis import (
    load_results_from_csv,
    filter_results_by_percentile,
    plot_score_distribution
)

def main():
    parser = argparse.ArgumentParser(description="Filtraggio dei risultati dell'analisi dei peptidi")
    
    parser.add_argument("--results-csv", required=True,
                        help="File CSV con i risultati dell'analisi")
    parser.add_argument("--output-dir", default="filtered_output",
                        help="Directory di output per i risultati filtrati")
    parser.add_argument("--percentile-threshold", type=float, default=1.0,
                        help="Soglia di percentile rank per il filtraggio")
    parser.add_argument("--percentile-operator", default="<", choices=["<", "<=", ">", ">=", "=="],
                        help="Operatore di confronto per il filtraggio del percentile rank")
    
    args = parser.parse_args()
    
    # Creiamo la directory di output se non esiste
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Carichiamo i risultati
    print(f"Caricamento dei risultati da {args.results_csv}...")
    results = load_results_from_csv(args.results_csv)
    
    if not results:
        print(f"Nessun risultato trovato in {args.results_csv}")
        return
    
    print(f"Caricati {len(results)} risultati.")
    
    # Filtriamo i risultati
    filtered_results = filter_results_by_percentile(
        results, args.percentile_threshold, args.percentile_operator
    )
    
    # Salviamo i risultati filtrati in un nuovo file CSV
    output_csv = os.path.join(args.output_dir, "filtered_results.csv")
    df = pd.DataFrame(filtered_results)
    # Ordiniamo per percentile_rank prima di salvare
    df = df.sort_values(by='percentile_rank')
    df.to_csv(output_csv, index=False)
    print(f"Risultati filtrati salvati in {output_csv} (ordinati per percentile_rank crescente)")
    
    # Creiamo un grafico della distribuzione degli score per i risultati filtrati
    output_plot = os.path.join(args.output_dir, "filtered_score_distribution.png")
    plot_score_distribution(filtered_results, output_plot)
    
    # Creiamo un istogramma dei percentile rank
    plt.figure(figsize=(10, 6))
    df = pd.DataFrame(filtered_results)
    plt.hist(df['percentile_rank'], bins=50, alpha=0.75)
    plt.title(f'Distribuzione dei percentile rank (filtrati con {args.percentile_operator} {args.percentile_threshold})')
    plt.xlabel('Percentile Rank')
    plt.ylabel('Frequenza')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    hist_output = os.path.join(args.output_dir, "percentile_rank_histogram.png")
    plt.savefig(hist_output, dpi=300, bbox_inches='tight')
    print(f"Istogramma dei percentile rank salvato in {hist_output}")
    
    # Stampiamo alcune statistiche sui risultati filtrati
    if filtered_results:
        df = pd.DataFrame(filtered_results)
        print("\nStatistiche sui risultati filtrati:")
        print(f"- Numero di peptidi: {len(df['peptide'].unique())}")
        print(f"- Numero di alleli: {len(df['allele'].unique())}")
        print(f"- Percentile rank minimo: {df['percentile_rank'].min():.4f}")
        print(f"- Percentile rank massimo: {df['percentile_rank'].max():.4f}")
        print(f"- Percentile rank medio: {df['percentile_rank'].mean():.4f}")
        print(f"- Score minimo: {df['score'].min():.4f}")
        print(f"- Score massimo: {df['score'].max():.4f}")
        print(f"- Score medio: {df['score'].mean():.4f}")

if __name__ == "__main__":
    main()
