#!/usr/bin/env python3
"""
Esempio di utilizzo della funzionalità di predizione dell'immunogenicità.
"""

import os
import argparse
import sys

# Aggiungiamo la directory principale al path per importare il pacchetto
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from peptide_analysis import (
    predict_peptide_immunogenicity,
    get_available_immunogenicity_alleles,
    print_available_immunogenicity_alleles,
    analyze_peptide_immunogenicity,
    load_peptides_from_csv
)

def main():
    parser = argparse.ArgumentParser(description="Predizione dell'immunogenicità dei peptidi")
    
    parser.add_argument("--input-csv", help="File CSV contenente i peptidi da analizzare")
    parser.add_argument("--output-dir", default="immunogenicity_output",
                        help="Directory di output per i risultati")
    parser.add_argument("--custom-mask", help="Maschera personalizzata (es. '2,3,9')")
    parser.add_argument("--allele", help="Allele HLA da utilizzare (es. 'HLA-A0201')")
    parser.add_argument("--list-alleles", action="store_true",
                        help="Mostra la lista degli alleli disponibili")
    parser.add_argument("--peptides", nargs="+", help="Lista di peptidi da analizzare")
    
    args = parser.parse_args()
    
    # Mostriamo la lista degli alleli disponibili se richiesto
    if args.list_alleles:
        print_available_immunogenicity_alleles()
        return
    
    # Verifichiamo che sia stato specificato almeno un input
    if not args.input_csv and not args.peptides:
        print("Errore: specificare un file CSV di input o una lista di peptidi.")
        parser.print_help()
        return
    
    # Creiamo la directory di output se non esiste
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Determiniamo i peptidi da analizzare
    peptides = None
    if args.peptides:
        peptides = args.peptides
    
    # Eseguiamo l'analisi
    report = analyze_peptide_immunogenicity(
        peptides=peptides,
        input_csv=args.input_csv,
        custom_mask=args.custom_mask,
        allele=args.allele,
        output_dir=args.output_dir
    )
    
    # Stampiamo un riepilogo dei risultati
    print("\nRiepilogo dell'analisi:")
    for key, value in report.items():
        if key != "error":
            print(f"- {key}: {value}")

if __name__ == "__main__":
    main()
