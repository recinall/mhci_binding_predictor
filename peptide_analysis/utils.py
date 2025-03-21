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

def filter_results_by_percentile(results, threshold=0.5, operator="<"):
    """
    Filtra i risultati in base al percentile rank.
    
    Parametri:
    results (list): Lista di risultati
    threshold (float): Soglia di percentile rank (default: 0.5)
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
        print(f"Operatore '{operator}' non valido. Utilizzo dell'operatore '<'")
        filtered = [r for r in results if r['percentile_rank'] < threshold]
        print(f"Filtrati {len(filtered)} risultati con percentile rank < {threshold}")
    
    return filtered

def filter_results_by_immunogenicity(results, threshold=0.0, operator=">"):
    """
    Filtra i risultati in base allo score di immunogenicità.
    
    Parametri:
    results (list): Lista di risultati
    threshold (float): Soglia di immunogenicità (default: 0.0)
    operator (str): Operatore di confronto ("<", "<=", ">", ">=", "==")
    
    Returns:
    list: Lista di risultati filtrati
    """
    # Verifichiamo che i risultati contengano lo score di immunogenicità
    if not results or 'immunogenicity_score' not in results[0]:
        print("I risultati non contengono lo score di immunogenicità.")
        return results
    
    # Filtriamo solo i risultati che hanno un valore di immunogenicità (non None)
    valid_results = [r for r in results if r['immunogenicity_score'] is not None]
    
    if operator == "<":
        filtered = [r for r in valid_results if r['immunogenicity_score'] < threshold]
        print(f"Filtrati {len(filtered)} risultati con immunogenicity score < {threshold}")
    elif operator == "<=":
        filtered = [r for r in valid_results if r['immunogenicity_score'] <= threshold]
        print(f"Filtrati {len(filtered)} risultati con immunogenicity score <= {threshold}")
    elif operator == ">":
        filtered = [r for r in valid_results if r['immunogenicity_score'] > threshold]
        print(f"Filtrati {len(filtered)} risultati con immunogenicity score > {threshold}")
    elif operator == ">=":
        filtered = [r for r in valid_results if r['immunogenicity_score'] >= threshold]
        print(f"Filtrati {len(filtered)} risultati con immunogenicity score >= {threshold}")
    elif operator == "==":
        filtered = [r for r in valid_results if r['immunogenicity_score'] == threshold]
        print(f"Filtrati {len(filtered)} risultati con immunogenicity score == {threshold}")
    else:
        print(f"Operatore '{operator}' non valido. Utilizzo dell'operatore '>'")
        filtered = [r for r in valid_results if r['immunogenicity_score'] > threshold]
        print(f"Filtrati {len(filtered)} risultati con immunogenicity score > {threshold}")
    
    return filtered

def filter_results_by_ic50(results, threshold=500, operator="<"):
    """
    Filtra i risultati in base al valore IC50.
    
    Parametri:
    results (list): Lista di risultati
    threshold (float): Soglia di IC50 (default: 500)
    operator (str): Operatore di confronto ("<", "<=", ">", ">=", "==")
    
    Returns:
    list: Lista di risultati filtrati
    """
    # Verifichiamo che i risultati contengano il valore IC50
    if not results:
        return results
    
    # Cerchiamo il campo IC50 nei risultati
    ic50_field = None
    for field in results[0].keys():
        if 'ic50' in field.lower():
            ic50_field = field
            break
    
    if not ic50_field:
        print("I risultati non contengono il valore IC50.")
        return results
    
    # Filtriamo solo i risultati che hanno un valore di IC50 (non None)
    valid_results = [r for r in results if r[ic50_field] is not None]
    
    if operator == "<":
        filtered = [r for r in valid_results if r[ic50_field] < threshold]
        print(f"Filtrati {len(filtered)} risultati con IC50 < {threshold}")
    elif operator == "<=":
        filtered = [r for r in valid_results if r[ic50_field] <= threshold]
        print(f"Filtrati {len(filtered)} risultati con IC50 <= {threshold}")
    elif operator == ">":
        filtered = [r for r in valid_results if r[ic50_field] > threshold]
        print(f"Filtrati {len(filtered)} risultati con IC50 > {threshold}")
    elif operator == ">=":
        filtered = [r for r in valid_results if r[ic50_field] >= threshold]
        print(f"Filtrati {len(filtered)} risultati con IC50 >= {threshold}")
    elif operator == "==":
        filtered = [r for r in valid_results if r[ic50_field] == threshold]
        print(f"Filtrati {len(filtered)} risultati con IC50 == {threshold}")
    else:
        print(f"Operatore '{operator}' non valido. Utilizzo dell'operatore '<'")
        filtered = [r for r in valid_results if r[ic50_field] < threshold]
        print(f"Filtrati {len(filtered)} risultati con IC50 < {threshold}")
    
    return filtered

def filter_results_combined(results, percentile_threshold=0.5, percentile_operator="<", 
                           immunogenicity_threshold=0.0, immunogenicity_operator=">",
                           ic50_threshold=500, ic50_operator="<"):
    """
    Filtra i risultati in base a percentile rank, score di immunogenicità e IC50.
    
    Parametri:
    results (list): Lista di risultati
    percentile_threshold (float): Soglia di percentile rank (default: 0.5)
    percentile_operator (str): Operatore di confronto per percentile ("<", "<=", ">", ">=", "==")
    immunogenicity_threshold (float): Soglia di immunogenicità (default: 0.0)
    immunogenicity_operator (str): Operatore di confronto per immunogenicità ("<", "<=", ">", ">=", "==")
    ic50_threshold (float): Soglia di IC50 (default: 500)
    ic50_operator (str): Operatore di confronto per IC50 ("<", "<=", ">", ">=", "==")
    
    Returns:
    list: Lista di risultati filtrati
    """
    # Prima filtriamo per percentile rank
    filtered_by_percentile = filter_results_by_percentile(
        results, percentile_threshold, percentile_operator
    )
    
    # Poi filtriamo per immunogenicità
    filtered_by_immunogenicity = filter_results_by_immunogenicity(
        filtered_by_percentile, immunogenicity_threshold, immunogenicity_operator
    )
    
    # Infine filtriamo per IC50
    filtered_combined = filter_results_by_ic50(
        filtered_by_immunogenicity, ic50_threshold, ic50_operator
    )
    
    print(f"Filtrati {len(filtered_combined)} risultati con tutti i criteri:")
    print(f"- Percentile rank {percentile_operator} {percentile_threshold}")
    print(f"- Immunogenicity score {immunogenicity_operator} {immunogenicity_threshold}")
    print(f"- IC50 {ic50_operator} {ic50_threshold}")
    
    return filtered_combined

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
