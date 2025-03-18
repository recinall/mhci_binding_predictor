#!/usr/bin/env python3
"""
Esempio di utilizzo della funzionalità di generazione e analisi di varianti di sequenze.
"""

import os
import argparse
import sys

# Aggiungiamo la directory principale al path per importare il pacchetto
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from peptide_analysis import (
    generate_sequence_variants,
    save_variants_to_csv,
    generate_all_variants,
    analyze_sequence_variants
)

# Sequenze di esempio (prese da pyIEDB.py)
DEFAULT_SEQUENCES = [
    [["S"], ["D","G","N"], ["P"], ["A","K","V"], ["R","C"], ["Y","A","H","N"], ["E","P","H"], ["F","H"], ["L"]],
    [["F","H"], ["L"], ["W","I","N"], ["G","V","L"], ["P","T","H"], ["R","N"], ["A","T"], ["L","V","H"], ["A","V","I"]],
    [["I"], ["F"], ["S","G"], ["K","T","E"], ["A","I"], ["S","A"], ["E","S","K"], ["S","Y","C"], ["L","M","F"]],
    [["K"], ["V","H"], ["A","K","E"], ["E","K","A"], ["L","S"], ["V","A","E"], ["H","K"], ["F","I"], ["L","F"]]
]

def main():
    parser = argparse.ArgumentParser(description="Generazione e analisi di varianti di sequenze peptidiche")
    
    parser.add_argument("--mode", choices=["generate", "analyze", "both"], default="both",
                        help="Modalità di esecuzione")
    parser.add_argument("--output-dir", default="sequence_variants_output",
                        help="Directory di output")
    parser.add_argument("--alleles", default="HLA-A*01:01,HLA-A*02:01",
                        help="Lista di alleli HLA separati da virgola")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Dimensione del batch per le richieste")
    parser.add_argument("--percentile-threshold", type=float, default=10.0,
                        help="Soglia di percentile rank per il filtraggio")
    parser.add_argument("--percentile-operator", default="<=", choices=["<", "<=", ">", ">=", "=="],
                        help="Operatore di confronto per il filtraggio del percentile rank")
    
    args = parser.parse_args()
    
    # Creiamo la directory di output se non esiste
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.mode == "generate" or args.mode == "both":
        print("Generazione delle varianti di sequenze...")
        total_variants, csv_files = generate_all_variants(DEFAULT_SEQUENCES, args.output_dir)
        print(f"Generate {total_variants} varianti in totale.")
        print(f"File CSV generati: {csv_files}")
    
    if args.mode == "analyze" or args.mode == "both":
        print("\nAnalisi delle varianti di sequenze...")
        report = analyze_sequence_variants(
            DEFAULT_SEQUENCES,
            args.alleles,
            args.batch_size,
            args.output_dir,
            args.percentile_threshold,
            args.percentile_operator
        )
        
        print("\nReport dell'analisi:")
        for key, value in report.items():
            if key != "csv_files":  # Evitiamo di stampare la lista completa dei file
                print(f"- {key}: {value}")

if __name__ == "__main__":
    main()
