#!/usr/bin/env python3
"""
Esempio di utilizzo della funzionalità di filtraggio combinato dei risultati.
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
    filter_results_by_immunogenicity,
    filter_results_combined,
    plot_score_distribution,
    plot_immunogenicity_correlation
)

def main():
    parser = argparse.ArgumentParser(description="Filtraggio combinato dei risultati dell'analisi dei peptidi")
    
    parser.add_argument("--results-csv", required=True,
                        help="File CSV con i risultati dell'analisi")
    parser.add_argument("--output-dir", default="filtered_combined_output",
                        help="Directory di output per i risultati filtrati")
    parser.add_argument("--percentile-threshold", type=float, default=0.5,
                        help="Soglia di percentile rank per il filtraggio")
    parser.add_argument("--percentile-operator", default="<", choices=["<", "<=", ">", ">=", "=="],
                        help="Operatore di confronto per il filtraggio del percentile rank")
    parser.add_argument("--immunogenicity-threshold", type=float, default=0.0,
                        help="Soglia di immunogenicità per il filtraggio")
    parser.add_argument("--immunogenicity-operator", default=">", choices=["<", "<=", ">", ">=", "=="],
                        help="Operatore di confronto per il filtraggio dell'immunogenicità")
    
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
    
    # Verifichiamo che i risultati contengano lo score di immunogenicità
    if 'immunogenicity_score' not in results[0]:
        print("I risultati non contengono lo score di immunogenicità.")
        print("Utilizzare prima run_complete_analysis con include_immunogenicity=True.")
        return
    
    # Filtriamo i risultati con entrambi i criteri
    filtered_results = filter_results_combined(
        results, 
        args.percentile_threshold, 
        args.percentile_operator,
        args.immunogenicity_threshold,
        args.immunogenicity_operator
    )
    
    # Salviamo i risultati filtrati in un nuovo file CSV
    output_csv = os.path.join(args.output_dir, "filtered_combined_results.csv")
    df = pd.DataFrame(filtered_results)
    # Ordiniamo per punteggio_composito decrescente
    if 'punteggio_composito' in df.columns:
        df = df.sort_values(by='punteggio_composito', ascending=False)
    else:
        # Altrimenti ordiniamo per percentile_rank crescente
        df = df.sort_values(by='percentile_rank')
    df.to_csv(output_csv, index=False)
    print(f"Risultati filtrati salvati in {output_csv}")
    
    # Creiamo un grafico della distribuzione degli score per i risultati filtrati
    output_plot = os.path.join(args.output_dir, "filtered_score_distribution.png")
    plot_score_distribution(filtered_results, output_plot)
    
    # Creiamo un grafico di correlazione tra percentile rank e immunogenicità
    output_corr = os.path.join(args.output_dir, "immunogenicity_correlation.png")
    plot_immunogenicity_correlation(filtered_results, output_corr)
    
    # Creiamo un grafico a dispersione di percentile rank vs immunogenicità
    plt.figure(figsize=(10, 8))
    df = pd.DataFrame(filtered_results)
    plt.scatter(df['percentile_rank'], df['immunogenicity_score'], alpha=0.7)
    plt.title('Percentile Rank vs Immunogenicità (risultati filtrati)')
    plt.xlabel('Percentile Rank')
    plt.ylabel('Immunogenicity Score')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    scatter_output = os.path.join(args.output_dir, "percentile_vs_immunogenicity.png")
    plt.savefig(scatter_output, dpi=300, bbox_inches='tight')
    print(f"Grafico a dispersione salvato in {scatter_output}")
    
    # Stampiamo alcune statistiche sui risultati filtrati
    if filtered_results:
        df = pd.DataFrame(filtered_results)
        print("\nStatistiche sui risultati filtrati:")
        print(f"- Numero di peptidi: {len(df['peptide'].unique())}")
        print(f"- Numero di alleli: {len(df['allele'].unique())}")
        print(f"- Percentile rank minimo: {df['percentile_rank'].min():.4f}")
        print(f"- Percentile rank massimo: {df['percentile_rank'].max():.4f}")
        print(f"- Percentile rank medio: {df['percentile_rank'].mean():.4f}")
        print(f"- Immunogenicity score minimo: {df['immunogenicity_score'].min():.4f}")
        print(f"- Immunogenicity score massimo: {df['immunogenicity_score'].max():.4f}")
        print(f"- Immunogenicity score medio: {df['immunogenicity_score'].mean():.4f}")
        
        # Stampiamo la distribuzione delle categorie se presenti
        if 'categoria' in df.columns:
            category_counts = df['categoria'].value_counts()
            print("\nDistribuzione delle categorie:")
            for category, count in category_counts.items():
                print(f"- {category}: {count}")

if __name__ == "__main__":
    main()
