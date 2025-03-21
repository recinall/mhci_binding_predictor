#!/usr/bin/env python
"""
Esempio di utilizzo del modulo binding_prediction per predire il binding MHC-I.
"""

import os
import sys
import pandas as pd

# Aggiungi la directory principale al path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from peptide_analysis.binding_prediction import predict_binding, get_available_methods, get_available_alleles

def main():
    """Funzione principale per l'esempio di predizione del binding."""
    
    print("Esempio di predizione del binding MHC-I")
    print("======================================")
    
    # Esempio di peptidi
    peptides = [
        "YLQPRTFLL",  # Peptide da SARS-CoV-2
        "KLPDDFTGCV", # Peptide da SARS-CoV-2
        "CINGVCWTV",  # Peptide da SARS-CoV-2
        "GILGFVFTL",  # Peptide da influenza
        "NLVPMVATV"   # Peptide da CMV
    ]
    
    # Filtra i peptidi per lunghezza (devono essere tutti della stessa lunghezza)
    peptides_9mer = [p for p in peptides if len(p) == 9]
    
    print(f"\nPredizione per {len(peptides_9mer)} peptidi di 9 amminoacidi:")
    for p in peptides_9mer:
        print(f"- {p}")
    
    # Predizione con netmhcpan_el (IEDB_recommended_epitope)
    try:
        print("\nPredizione con netmhcpan_el (IEDB_recommended_epitope)...")
        results = predict_binding(
            peptides=peptides_9mer,
            allele="HLA-A*02:01",
            method="netmhcpan_el"
        )
        
        # Mostra i risultati
        if not results.empty:
            print("\nRisultati della predizione:")
            pd.set_option('display.max_columns', None)
            print(results)
            
            # Ordina per rank (migliore in cima)
            if 'rank' in results.columns:
                sorted_results = results.sort_values(by='rank')
                print("\nRisultati ordinati per rank (migliore in cima):")
                print(sorted_results[['peptide', 'rank']])
        else:
            print("Nessun risultato ottenuto.")
    
    except Exception as e:
        print(f"Errore durante la predizione: {e}")
        
        # Prova a ottenere i metodi disponibili
        try:
            print("\nMetodi di predizione disponibili:")
            methods = get_available_methods()
            for method in methods:
                print(f"- {method}")
        except Exception as e2:
            print(f"Impossibile ottenere i metodi disponibili: {e2}")

if __name__ == "__main__":
    main()
