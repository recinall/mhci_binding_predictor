"""
Modulo di utilità per la libreria di analisi dei peptidi.
"""

import csv
import pandas as pd
import os

def load_peptides_from_csv(csv_file):
    """
    Carica i peptidi da un file CSV.
    
    Parametri:
    csv_file (str): Nome del file CSV
    
    Returns:
    list: Lista di peptidi
    """
    peptides = []
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                peptides.append(row['peptide'])
        
        print(f"Caricati {len(peptides)} peptidi da {csv_file}")
        return peptides
    except Exception as e:
        print(f"Errore nel caricamento dei peptidi da {csv_file}: {e}")
        return []

def load_results_from_csv(csv_file, sort_by_percentile=True):
    """
    Carica i risultati dell'analisi da un file CSV.
    
    Parametri:
    csv_file (str): Nome del file CSV
    sort_by_percentile (bool): Se True, ordina i risultati per percentile_rank crescente
    
    Returns:
    list: Lista di risultati
    """
    results = []
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Convertiamo i valori numerici
                result = {
                    'peptide': row['peptide'],
                    'allele': row['allele'],
                    'score': float(row['score']),
                    'percentile_rank': float(row['percentile_rank'])
                }
                results.append(result)
        
        # Ordiniamo i risultati per percentile_rank se richiesto
        if sort_by_percentile and results:
            results = sorted(results, key=lambda x: x['percentile_rank'])
            print(f"Caricati e ordinati {len(results)} risultati da {csv_file}")
        else:
            print(f"Caricati {len(results)} risultati da {csv_file}")
        
        return results
    except Exception as e:
        print(f"Errore nel caricamento dei risultati da {csv_file}: {e}")
        return []

def ensure_directory(directory):
    """
    Assicura che una directory esista, creandola se necessario.
    
    Parametri:
    directory (str): Percorso della directory
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Creata directory: {directory}")

def filter_results_by_percentile(results, threshold=1.0, operator="<="):
    """
    Filtra i risultati in base al percentile rank.
    
    Parametri:
    results (list): Lista di risultati
    threshold (float): Soglia di percentile rank (default: 1.0)
    operator (str): Operatore di confronto ("<", "<=", ">", ">=", "==")
    
    Returns:
    list: Lista di risultati filtrati
    """
    if operator == "<":
        filtered = [r for r in results if r['percentile_rank'] < threshold]
        print(f"Filtrati {len(filtered)} risultati con percentile rank < {threshold}")
    elif operator == "<=":
        filtered = [r for r in results if r['percentile_rank'] <= threshold]
        print(f"Filtrati {len(filtered)} risultati con percentile rank <= {threshold}")
    elif operator == ">":
        filtered = [r for r in results if r['percentile_rank'] > threshold]
        print(f"Filtrati {len(filtered)} risultati con percentile rank > {threshold}")
    elif operator == ">=":
        filtered = [r for r in results if r['percentile_rank'] >= threshold]
        print(f"Filtrati {len(filtered)} risultati con percentile rank >= {threshold}")
    elif operator == "==":
        filtered = [r for r in results if r['percentile_rank'] == threshold]
        print(f"Filtrati {len(filtered)} risultati con percentile rank == {threshold}")
    else:
        print(f"Operatore '{operator}' non valido. Utilizzo dell'operatore '<='")
        filtered = [r for r in results if r['percentile_rank'] <= threshold]
        print(f"Filtrati {len(filtered)} risultati con percentile rank <= {threshold}")
    
    return filtered

def export_results_to_excel(results, output_file):
    """
    Esporta i risultati in un file Excel.
    
    Parametri:
    results (list): Lista di risultati
    output_file (str): Nome del file Excel
    """
    df = pd.DataFrame(results)
    
    try:
        df.to_excel(output_file, index=False)
        print(f"Risultati esportati in {output_file}")
    except Exception as e:
        print(f"Errore nell'esportazione dei risultati in {output_file}: {e}")
