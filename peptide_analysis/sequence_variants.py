"""
Modulo per la generazione di varianti di sequenze peptidiche.
"""

import csv
import os
from tqdm import tqdm

def generate_sequence_variants(sequenza):
    """
    Genera tutte le permutazioni possibili di una sequenza di amminoacidi,
    dove alcune posizioni possono avere varianti.
    
    Parametri:
    sequenza (list): Una lista di liste, dove ogni lista interna contiene 
                    gli amminoacidi possibili per quella posizione.
                    Esempio: [['A'], ['B', 'C'], ['D', 'E', 'F']]
    
    Returns:
    list: Lista di tutte le possibili stringhe permutate
    """
    risultato = [""]
    
    for posizione in sequenza:
        nuove_combinazioni = []
        for stringa_parziale in risultato:
            for variante in posizione:
                nuove_combinazioni.append(stringa_parziale + variante)
        risultato = nuove_combinazioni
    
    return risultato

def save_variants_to_csv(peptides, filename="sequence_variants.csv"):
    """
    Salva i peptidi varianti in un file CSV.
    
    Parametri:
    peptides (list): Lista di peptidi
    filename (str): Nome del file CSV
    
    Returns:
    str: Nome del file CSV creato
    """
    print(f"Salvando {len(peptides)} varianti peptidiche in {filename}...")
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['peptide'])  # Header
        for peptide in peptides:
            writer.writerow([peptide])
    
    print(f"Salvato con successo in {filename}")
    return filename

def generate_all_variants(sequences, output_dir="variants", combined_file="Sequenze_AIO.csv"):
    """
    Genera tutte le varianti per un insieme di sequenze e le salva in file CSV separati.
    
    Parametri:
    sequences (list): Lista di sequenze, dove ogni sequenza è una lista di liste
    output_dir (str): Directory di output per i file CSV
    combined_file (str): Nome del file CSV combinato con tutte le varianti
    
    Returns:
    tuple: (numero totale di varianti, lista di file CSV generati)
    """
    # Assicuriamo che la directory di output esista
    os.makedirs(output_dir, exist_ok=True)
    
    all_variants = []
    csv_files = []
    
    print(f"Generazione di varianti per {len(sequences)} sequenze...")
    
    for i, sequence in enumerate(sequences):
        print(f"Sequenza {i+1}: generazione varianti...")
        variants = generate_sequence_variants(sequence)
        print(f"Sequenza {i+1}: generate {len(variants)} varianti")
        
        # Salviamo le varianti in un file CSV
        filename = os.path.join(output_dir, f"Sequenza_{i+1:02d}.csv")
        save_variants_to_csv(variants, filename)
        csv_files.append(filename)
        
        # Aggiungiamo le varianti alla lista completa
        all_variants.extend(variants)
    
    # Salviamo tutte le varianti in un unico file
    combined_path = os.path.join(output_dir, combined_file)
    save_variants_to_csv(all_variants, combined_path)
    csv_files.append(combined_path)
    
    print(f"Generazione completata. Totale varianti: {len(all_variants)}")
    print(f"File generati: {len(csv_files)}")
    
    return len(all_variants), csv_files
