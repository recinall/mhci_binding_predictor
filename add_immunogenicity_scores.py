#!/usr/bin/env python

import csv
import os,sys
import tempfile
from immunogenicity import predict_immunogenicity

Prediction = predict_immunogenicity.Prediction

def add_immunogenicity_scores(input_csv, output_csv=None, filtered_csv=None, ranked_csv=None):
    """
    Aggiunge il punteggio di immunogenicità calcolato con la classe Prediction
    al file CSV contenente peptidi e altri punteggi, considerando l'allele specifico.
    
    Args:
        input_csv: Path del file CSV di input
        output_csv: Path del file CSV di output (opzionale)
        filtered_csv: Path del file CSV filtrato (opzionale)
        ranked_csv: Path del file CSV con peptidi classificati (opzionale)
    """
    if output_csv is None:
        base_name = os.path.splitext(input_csv)[0]
        output_csv = f"{base_name}_with_immunogenicity.csv"
    
    if filtered_csv is None:
        base_name = os.path.splitext(input_csv)[0]
        filtered_csv = f"{base_name}_filtered.csv"
        
    if ranked_csv is None:
        base_name = os.path.splitext(input_csv)[0]
        ranked_csv = f"{base_name}_ranked.csv"
    
    # Legge il file CSV di input
    with open(input_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
    
    # Inizializza la classe Prediction
    pred = Prediction()
    
    # Raggruppa i peptidi per allele
    peptides_by_allele = {}
    for row in rows:
        allele = row['allele']
        peptide = row['peptide']
        if allele not in peptides_by_allele:
            peptides_by_allele[allele] = []
        peptides_by_allele[allele].append(peptide)
    
    # Calcola l'immunogenicità per ogni gruppo di peptidi con lo stesso allele
    immunogenicity_scores = {}
    
    for allele, peptides in peptides_by_allele.items():
        # Crea un file temporaneo contenente i peptidi (uno per riga)
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            for peptide in peptides:
                temp_file.write(f"{peptide}\n")
            temp_file_path = temp_file.name
        
        # Simula le opzioni necessarie per la funzione validate
        class Options:
            def __init__(self, allele_value):
                self.custom_mask = None
                self.allele = allele_value
        
        # Rimuovi simboli * e : dall'allele per corrispondere al formato richiesto dalla classe Prediction
        allele_formatted = allele.replace("*", "").replace(":", "")
        options = Options(allele_formatted)
        
        try:
            # Ottiene i dati validati
            peptides_validated, custom_mask, allele_used = pred.validate(options, [temp_file_path])
            
            # Calcola i punteggi di immunogenicità
            immunoscale = {"A":0.127, "C":-0.175, "D":0.072, "E":0.325, "F":0.380, "G":0.110, "H":0.105, "I":0.432, "K":-0.700, "L":-0.036, "M":-0.570, "N":-0.021, "P":-0.036, "Q":-0.376, "R":0.168, "S":-0.537, "T":0.126, "V":0.134, "W":0.719, "Y":-0.012}
            immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
            
            # Determina la mask in base all'allele
            if allele_used and custom_mask:
                mask_str = custom_mask.split(",")
                mask_num = list(map(int, mask_str))
                mask_num = list(map(lambda x: x - 1, mask_num))  # Convert to 0-based indexing
            else:
                # Default mask (1st, 2nd, and C-terminus)
                for peptide in peptides_validated:
                    peplen = len(peptide)
                    cterm = peplen - 1
                    mask_num = [0, 1, cterm]
            
            # Calcola il punteggio per ogni peptide
            for peptide in peptides_validated:
                peplen = len(peptide)
                score = 0
                count = 0
                
                # Adatta i pesi in base alla lunghezza del peptide
                if peplen > 9:
                    pepweight = immunoweight[:5] + ((peplen - 9) * [0.30]) + immunoweight[5:]
                else:
                    pepweight = immunoweight
                
                # Calcola il punteggio
                for pos in peptide:
                    if count not in mask_num:
                        score += pepweight[count] * immunoscale[pos]
                    count += 1
                
                immunogenicity_scores[(peptide, allele)] = round(score, 5)
        
        except Exception as e:
            print(f"Errore durante il calcolo dell'immunogenicità per l'allele {allele}: {str(e)}")
            # Assegna un valore predefinito o marca come errore
            for peptide in peptides:
                immunogenicity_scores[(peptide, allele)] = "Error"
        
        finally:
            # Rimuove il file temporaneo
            os.unlink(temp_file_path)
    
    # Aggiunge i punteggi di immunogenicità alle righe del CSV
    for row in rows:
        peptide = row['peptide']
        allele = row['allele']
        row['immunogenicity'] = immunogenicity_scores.get((peptide, allele), "N/A")
    
    # Scrive il file CSV di output completo
    fieldnames = list(rows[0].keys())
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"File CSV con punteggi di immunogenicità creato: {output_csv}")
    
    # Filtra le righe in base ai criteri: percentile_rank < 0.5 e immunogenicity > 0
    filtered_rows = []
    for row in rows:
        try:
            percentile_rank = float(row['percentile_rank'])
            immunogenicity = float(row['immunogenicity'])
            
            if percentile_rank < 0.5 and immunogenicity > 0:
                filtered_rows.append(row)
        except (ValueError, TypeError):
            # Salta le righe con valori non numerici
            continue
    
    # Scrive il file CSV filtrato
    if filtered_rows:
        with open(filtered_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(filtered_rows)
        
        print(f"File CSV filtrato creato: {filtered_csv}")
        print(f"Numero di peptidi nel file filtrato: {len(filtered_rows)}")
    else:
        print("Nessun peptide soddisfa i criteri di filtro.")
    
    # Classifica i peptidi
    ranked_peptides = rank_peptides(rows)
    
    # Scrive il file CSV con le classificazioni
    if ranked_peptides:
        with open(ranked_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['peptide', 'allele', 'composite_score', 'category', 'immunogenicity', 'percentile_rank', 'score'])
            writer.writerows(ranked_peptides)
        
        print(f"File CSV con peptidi classificati creato: {ranked_csv}")
    
    return output_csv, filtered_csv, ranked_csv

def rank_peptides(rows):
    """
    Classifica i peptidi in base a una valutazione multifattoriale.
    
    Args:
        rows: Lista di dizionari rappresentanti le righe del CSV
    
    Returns:
        Lista di tuple (peptide, allele, punteggio_composito, categoria, immunogenicity, percentile_rank, score)
    """
    ranked_peptides = []
    
    for row in rows:
        try:
            peptide = row['peptide']
            allele = row['allele']
            percentile_rank = float(row['percentile_rank'])
            score = float(row['score'])
            immunogenicity = float(row['immunogenicity'])
            
            # Escludiamo subito i peptidi con immunogenicità negativa
            if immunogenicity <= 0:
                continue
            
            # Calcoliamo un punteggio composito
            composite_score = (immunogenicity * 0.5) + ((1 - percentile_rank/100) * 0.3) + (score * 0.2)
            
            # Assegniamo una categoria
            if immunogenicity > 0.3 and percentile_rank < 0.1 and score > 0.95:
                category = "Eccellente"
            elif immunogenicity > 0 and percentile_rank < 0.5 and score > 0.9:
                category = "Buono"
            elif immunogenicity > 0 and percentile_rank < 1.0 and score > 0.8:
                category = "Da considerare"
            else:
                category = "Non prioritario"
            
            ranked_peptides.append((peptide, allele, composite_score, category, immunogenicity, percentile_rank, score))
            
        except (ValueError, TypeError):
            continue
    
    # Ordiniamo per punteggio composito decrescente
    ranked_peptides.sort(key=lambda x: x[2], reverse=True)
    
    return ranked_peptides

if __name__ == "__main__":
    input_csv = sys.argv[1]
    add_immunogenicity_scores(input_csv)
