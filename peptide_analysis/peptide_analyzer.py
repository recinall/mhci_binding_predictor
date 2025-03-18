"""
Modulo per l'analisi di peptidi 9-mer utilizzando il servizio IEDB.
"""

import requests
import csv
import time
import pandas as pd
from io import StringIO
import os
from tqdm import tqdm

def send_iedb_request(sequence_text, allele):
    """
    Invia una richiesta al server IEDB per l'analisi MHC-I.
    
    Parametri:
    sequence_text (str): Sequenza peptidica da analizzare
    allele (str): Alleli HLA da utilizzare per la predizione
    
    Returns:
    str: Risposta del server come testo
    """
    url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"

    length = ",".join(["9" for _ in range(len(allele.split(",")))])
    payload = {
        "method": "recommended",
        "sequence_text": sequence_text,
        "allele": allele,
        "length": length
    }
    
    max_retries = 3
    retry_count = 0
    
    while retry_count < max_retries:
        try:
            response = requests.post(url, data=payload)
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            print(f"Errore nella richiesta: {e}")
            retry_count += 1
            if retry_count < max_retries:
                print(f"Ritento {retry_count}/{max_retries} dopo 5 secondi...")
                time.sleep(5)
            else:
                print("Numero massimo di tentativi raggiunto.")
                return None

def parse_iedb_response(response_text):
    """
    Parsa la risposta del server IEDB.
    
    Parametri:
    response_text (str): Testo della risposta del server
    
    Returns:
    list: Lista di dizionari con i risultati
    """
    if not response_text:
        return []
    
    # Utilizziamo pandas per parsare il testo tabulare
    df = pd.read_csv(StringIO(response_text), sep='\t')
    
    # Estraiamo solo le colonne che ci interessano
    results = []
    for _, row in df.iterrows():
        results.append({
            'peptide': row['peptide'],
            'allele': row['allele'],
            'score': row['score'],
            'percentile_rank': row['percentile_rank']
        })
    
    return results

def analyze_peptides(peptides, allele_list="HLA-A*01:01,HLA-A*02:01", batch_size=10, output_csv=None):
    """
    Analizza una lista di peptidi utilizzando il servizio IEDB.
    
    Parametri:
    peptides (list): Lista di peptidi da analizzare
    allele_list (str): Lista di alleli HLA da utilizzare
    batch_size (int): Dimensione del batch per le richieste
    output_csv (str): Nome del file CSV di output (opzionale)
    
    Returns:
    list: Lista di risultati dell'analisi
    """
    all_results = []
    
    # Processiamo i peptidi in batch
    total_peptides = len(peptides)
    
    with tqdm(total=total_peptides) as pbar:
        for i in range(0, total_peptides, batch_size):
            batch = peptides[i:i+batch_size]
            
            # Creiamo una stringa con i peptidi del batch
            sequence_text = "\n".join(batch)
            
            # Inviamo la richiesta al server
            response_text = send_iedb_request(sequence_text, allele_list)
            
            # Parsifichiamo la risposta
            if response_text:
                results = parse_iedb_response(response_text)
                all_results.extend(results)
            
            # Aggiorniamo la barra di progresso
            pbar.update(len(batch))
            
            # Pausa per evitare di sovraccaricare il server
            time.sleep(1)
    
    # Se è stato specificato un file di output, salviamo i risultati
    if output_csv:
        # Ordiniamo i risultati per percentile_rank (dal minimo al massimo)
        sorted_results = sorted(all_results, key=lambda x: x['percentile_rank'])
        
        with open(output_csv, 'w', newline='') as csvfile:
            fieldnames = ['peptide', 'allele', 'score', 'percentile_rank']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in sorted_results:
                writer.writerow(result)
        
        print(f"Risultati salvati in {output_csv} (ordinati per percentile_rank crescente)")
    
    return all_results

def process_peptides_from_csv(input_csv, output_csv, allele_list="HLA-A*01:01,HLA-A*02:01", batch_size=10):
    """
    Processa un file CSV contenente peptidi e salva i risultati in un nuovo file CSV.
    
    Parametri:
    input_csv (str): Nome del file CSV di input
    output_csv (str): Nome del file CSV di output
    allele_list (str): Lista di alleli HLA da utilizzare
    batch_size (int): Dimensione del batch per le richieste
    
    Returns:
    list: Lista di risultati dell'analisi
    """
    # Leggiamo il file CSV di input
    peptides = []
    with open(input_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            peptides.append(row['peptide'])
    
    # Analizziamo i peptidi
    results = analyze_peptides(peptides, allele_list, batch_size, output_csv)
    
    return results
