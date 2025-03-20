#!/usr/bin/env python3
"""
Esempio di utilizzo della classe CombinedResult per la gestione dei risultati combinati.
"""

import os
import argparse
import sys
import pandas as pd

# Aggiungiamo la directory principale al path per importare il pacchetto
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from peptide_analysis import (
    load_results_from_csv,
    add_immunogenicity_to_results,
    CombinedResult
)

def main():
    parser = argparse.ArgumentParser(description="Gestione dei risultati combinati di binding MHC-I e immunogenicità")
    
    parser.add_argument("--results-csv", required=True,
                        help="File CSV con i risultati dell'analisi di binding MHC-I")
    parser.add_argument("--output-dir", default="combined_results_output",
                        help="Directory di output per i risultati combinati")
    parser.add_argument("--create-excel", action="store_true",
                        help="Crea un file Excel con tutti i risultati")
    
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
    
    # Aggiungiamo l'immunogenicità ai risultati
    print("\nAggiunta dell'immunogenicità ai risultati...")
    results_with_immuno, filtered_results, ranked_results = add_immunogenicity_to_results(
        results, args.output_dir
    )
    
    # Creiamo un file Excel con tutti i risultati se richiesto
    if args.create_excel:
        excel_file = os.path.join(args.output_dir, "combined_results.xlsx")
        
        with pd.ExcelWriter(excel_file) as writer:
            # Foglio con tutti i risultati
            pd.DataFrame(results_with_immuno).to_excel(writer, sheet_name="Tutti i risultati", index=False)
            
            # Foglio con i risultati filtrati
            pd.DataFrame(filtered_results).to_excel(writer, sheet_name="Risultati filtrati", index=False)
            
            # Foglio con i risultati ordinati
            pd.DataFrame(ranked_results).to_excel(writer, sheet_name="Risultati ordinati", index=False)
            
            # Foglio con statistiche
            stats = {
                "Metrica": [
                    "Totale risultati",
                    "Risultati con immunogenicità",
                    "Risultati filtrati",
                    "Risultati ordinati",
                    "Percentile rank medio",
                    "Immunogenicità media",
                    "Punteggio composito medio"
                ],
                "Valore": [
                    len(results),
                    len(results_with_immuno),
                    len(filtered_results),
                    len(ranked_results),
                    round(pd.DataFrame(results_with_immuno)["percentile_rank"].mean(), 4),
                    round(pd.DataFrame(results_with_immuno)["immunogenicity_score"].mean(), 4),
                    round(pd.DataFrame(results_with_immuno)["punteggio_composito"].mean(), 4)
                ]
            }
            pd.DataFrame(stats).to_excel(writer, sheet_name="Statistiche", index=False)
        
        print(f"\nFile Excel con tutti i risultati creato: {excel_file}")
    
    # Stampiamo alcune statistiche sui risultati
    print("\nStatistiche sui risultati:")
    print(f"- Totale risultati: {len(results)}")
    print(f"- Risultati con immunogenicità: {len(results_with_immuno)}")
    print(f"- Risultati filtrati (percentile < 0.5 e immunogenicità > 0): {len(filtered_results)}")
    print(f"- Risultati ordinati per punteggio composito: {len(ranked_results)}")
    
    # Stampiamo la distribuzione delle categorie
    if results_with_immuno:
        df = pd.DataFrame(results_with_immuno)
        category_counts = df['categoria'].value_counts()
        print("\nDistribuzione delle categorie:")
        for category, count in category_counts.items():
            print(f"- {category}: {count}")

if __name__ == "__main__":
    main()
