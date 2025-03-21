#!/usr/bin/env python
"""
Esempio di utilizzo delle funzioni di filtraggio con IC50.
"""

import os
import sys
import pandas as pd

# Aggiungi la directory principale al path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from peptide_analysis.binding_prediction import predict_binding
from peptide_analysis.immunogenicity import predict_peptide_immunogenicity
from peptide_analysis.utils import (
    filter_results_by_percentile,
    filter_results_by_immunogenicity,
    filter_results_by_ic50,
    filter_results_combined
)

def main():
    """Funzione principale per l'esempio di filtraggio con IC50."""
    
    print("Esempio di filtraggio con IC50")
    print("=============================")
    
    # Esempio di peptidi
    peptides = [
        "YLQPRTFLL",  # Peptide da SARS-CoV-2
        "GILGFVFTL",  # Peptide da influenza
        "NLVPMVATV",  # Peptide da CMV
        "CINGVCWTV",  # Peptide da SARS-CoV-2
        "KLQPRTFLL"   # Variante di SARS-CoV-2
    ]
    
    # Predizione del binding
    try:
        print("\nPredizione del binding MHC-I...")
        binding_results = predict_binding(
            peptides=peptides,
            allele="HLA-A*02:01",
            method="netmhcpan_el"
        )
        
        if binding_results.empty:
            print("Nessun risultato di binding ottenuto.")
            return
        
        print("\nRisultati della predizione di binding:")
        print(binding_results[['peptide', 'score', 'rank']])
        
        # Predizione dell'immunogenicità
        print("\nPredizione dell'immunogenicità...")
        immunogenicity_results = predict_peptide_immunogenicity(
            peptides=peptides,
            allele="HLA-A*02:01"
        )
        
        # Combina i risultati
        combined_results = []
        for i, row in binding_results.iterrows():
            peptide = row['peptide']
            # Trova il risultato di immunogenicità corrispondente
            immuno_score = None
            for ir in immunogenicity_results:
                if ir['peptide'] == peptide:
                    immuno_score = ir['score']
                    break
            
            # Crea un dizionario con tutti i dati
            result = {
                'peptide': peptide,
                'allele': 'HLA-A*02:01',
                'percentile_rank': row['rank'],
                'score': row['score']
            }
            
            # Aggiungi IC50 se disponibile
            if 'ic50' in row:
                result['ic50'] = row['ic50']
            elif 'ann_ic50' in row:
                result['ic50'] = row['ann_ic50']
            
            # Aggiungi immunogenicity score
            if immuno_score is not None:
                result['immunogenicity_score'] = immuno_score
            
            combined_results.append(result)
        
        # Mostra i risultati combinati
        print("\nRisultati combinati:")
        for r in combined_results:
            print(f"Peptide: {r['peptide']}, Rank: {r['percentile_rank']:.4f}, "
                  f"Immunogenicità: {r.get('immunogenicity_score', 'N/A')}, "
                  f"IC50: {r.get('ic50', 'N/A')}")
        
        # Filtra per percentile rank
        print("\nFiltraggio per percentile rank < 0.5:")
        filtered_by_percentile = filter_results_by_percentile(
            combined_results, threshold=0.5, operator="<"
        )
        for r in filtered_by_percentile:
            print(f"Peptide: {r['peptide']}, Rank: {r['percentile_rank']:.4f}")
        
        # Filtra per immunogenicità
        print("\nFiltraggio per immunogenicità > 0:")
        filtered_by_immunogenicity = filter_results_by_immunogenicity(
            combined_results, threshold=0, operator=">"
        )
        for r in filtered_by_immunogenicity:
            print(f"Peptide: {r['peptide']}, Immunogenicità: {r.get('immunogenicity_score', 'N/A')}")
        
        # Filtra per IC50
        print("\nFiltraggio per IC50 < 500:")
        filtered_by_ic50 = filter_results_by_ic50(
            combined_results, threshold=500, operator="<"
        )
        for r in filtered_by_ic50:
            print(f"Peptide: {r['peptide']}, IC50: {r.get('ic50', 'N/A')}")
        
        # Filtra con tutti i criteri
        print("\nFiltraggio con tutti i criteri:")
        filtered_combined = filter_results_combined(
            combined_results,
            percentile_threshold=0.5, percentile_operator="<",
            immunogenicity_threshold=0, immunogenicity_operator=">",
            ic50_threshold=500, ic50_operator="<"
        )
        for r in filtered_combined:
            print(f"Peptide: {r['peptide']}, Rank: {r['percentile_rank']:.4f}, "
                  f"Immunogenicità: {r.get('immunogenicity_score', 'N/A')}, "
                  f"IC50: {r.get('ic50', 'N/A')}")
        
    except Exception as e:
        print(f"Errore durante l'esecuzione: {e}")

if __name__ == "__main__":
    main()
