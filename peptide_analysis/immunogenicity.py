"""
Modulo per la predizione dell'immunogenicità dei peptidi.

Questo modulo è basato sul codice originale di Dorjee Gyaltsen,
adattato per l'integrazione nella libreria peptide_analysis.
"""

import os
import sys
import csv
from tqdm import tqdm

class Immunogenicity:
    """
    Classe per la predizione dell'immunogenicità dei peptidi.
    
    Questa classe implementa l'algoritmo di predizione dell'immunogenicità
    descritto in Calis et al. (2013).
    """
    
    def __init__(self):
        """Inizializza la classe Immunogenicity."""
        # Dizionario degli alleli disponibili con le posizioni mascherate
        self.allele_dict = {
            "H-2-Db": "2,5,9", "H-2-Dd": "2,3,5", "H-2-Kb": "2,3,9", 
            "H-2-Kd": "2,5,9", "H-2-Kk": "2,8,9", "H-2-Ld": "2,5,9", 
            "HLA-A0101": "2,3,9", "HLA-A0201": "1,2,9", "HLA-A0202": "1,2,9", 
            "HLA-A0203": "1,2,9", "HLA-A0206": "1,2,9", "HLA-A0211": "1,2,9", 
            "HLA-A0301": "1,2,9", "HLA-A1101": "1,2,9", "HLA-A2301": "2,7,9", 
            "HLA-A2402": "2,7,9", "HLA-A2601": "1,2,9", "HLA-A2902": "2,7,9", 
            "HLA-A3001": "1,3,9", "HLA-A3002": "2,7,9", "HLA-A3101": "1,2,9", 
            "HLA-A3201": "1,2,9", "HLA-A3301": "1,2,9", "HLA-A6801": "1,2,9", 
            "HLA-A6802": "1,2,9", "HLA-A6901": "1,2,9", "HLA-B0702": "1,2,9", 
            "HLA-B0801": "2,5,9", "HLA-B1501": "1,2,9", "HLA-B1502": "1,2,9", 
            "HLA-B1801": "1,2,9", "HLA-B2705": "2,3,9", "HLA-B3501": "1,2,9", 
            "HLA-B3901": "1,2,9", "HLA-B4001": "1,2,9", "HLA-B4002": "1,2,9", 
            "HLA-B4402": "2,3,9", "HLA-B4403": "2,3,9", "HLA-B4501": "1,2,9", 
            "HLA-B4601": "1,2,9", "HLA-B5101": "1,2,9", "HLA-B5301": "1,2,9", 
            "HLA-B5401": "1,2,9", "HLA-B5701": "1,2,9", "HLA-B5801": "1,2,9"
        }
        
        # Scala di immunogenicità per gli amminoacidi
        self.immunoscale = {
            "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, 
            "G": 0.110, "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036, 
            "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168, 
            "S": -0.537, "T": 0.126, "V": 0.134, "W": 0.719, "Y": -0.012
        }
        
        # Pesi per le posizioni nel peptide
        self.immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
    
    def isint(self, x):
        """
        Verifica se un valore può essere convertito in un intero.
        
        Parametri:
        x: Valore da verificare
        
        Returns:
        bool: True se il valore può essere convertito in un intero, False altrimenti
        """
        try:
            a = float(x)
            b = int(a)
        except ValueError:
            return False
        else:
            return a == b
    
    def validate_peptides(self, peptides, custom_mask=None, allele=None):
        """
        Valida i peptidi e le opzioni di predizione.
        
        Parametri:
        peptides (list): Lista di peptidi da analizzare
        custom_mask (str): Maschera personalizzata (opzionale)
        allele (str): Allele HLA da utilizzare (opzionale)
        
        Returns:
        tuple: (peptidi validati, maschera personalizzata, allele)
        """
        # Validazione dei peptidi
        valid_peptides = []
        for peptide in peptides:
            peptide = peptide.strip().upper()
            is_valid = True
            
            for i, amino_acid in enumerate(peptide):
                if amino_acid not in "ACDEFGHIKLMNPQRSTVWY":
                    print(f"Peptide: '{peptide}' contiene un carattere non valido: '{amino_acid}' in posizione {i+1}.")
                    is_valid = False
                    break
            
            if is_valid:
                valid_peptides.append(peptide)
        
        # Validazione della maschera personalizzata
        if custom_mask:
            try:
                custom_mask_list = list(map(int, custom_mask.split(",")))
                if sum(n < 0 for n in custom_mask_list) > 0:
                    print("La maschera personalizzata deve contenere valori maggiori di zero.")
                    custom_mask = None
                
                max_length = max(custom_mask_list)
                if not all([len(peptide) >= max_length for peptide in valid_peptides]):
                    print(f"La lunghezza della maschera '{max_length}' non può essere maggiore della lunghezza del peptide.")
                    custom_mask = None
                
                if not all(self.isint(num) for num in custom_mask.split(",")):
                    print("La maschera personalizzata deve essere un numero singolo o una lista di numeri separati da virgola.")
                    custom_mask = None
            except:
                print("Errore nella validazione della maschera personalizzata.")
                custom_mask = None
        
        # Validazione dell'allele
        if allele:
            allele = allele.replace("*", "").replace(":", "")
            
            if allele in self.allele_dict:
                # Se è specificato sia l'allele che la maschera personalizzata, l'allele ha la precedenza
                if custom_mask:
                    print(f"* L'allele {allele} ha valore predefinito {self.allele_dict[allele]}.")
                    print(f"* Quando vengono utilizzate entrambe le opzioni 'custom_mask' e 'allele', quest'ultima ha la precedenza.")
                
                custom_mask = self.allele_dict[allele]
            else:
                print(f"L'allele {allele} non è disponibile.")
                allele = None
        
        return valid_peptides, custom_mask, allele
    
    def predict_immunogenicity(self, peptides, custom_mask=None, allele=None):
        """
        Predice l'immunogenicità per una lista di peptidi.
        
        Parametri:
        peptides (list): Lista di peptidi da analizzare
        custom_mask (str): Maschera personalizzata (opzionale)
        allele (str): Allele HLA da utilizzare (opzionale)
        
        Returns:
        list: Lista di risultati della predizione
        """
        # Validazione dei dati di input
        valid_peptides, custom_mask, allele = self.validate_peptides(peptides, custom_mask, allele)
        
        if not valid_peptides:
            print("Nessun peptide valido da analizzare.")
            return []
        
        result_list = []
        
        # Determiniamo la maschera da utilizzare
        if not custom_mask:
            mask_out = [1, 2, "cterm"]
        else:
            try:
                mask_str = custom_mask.split(",")
                mask_num = list(map(int, mask_str))
                mask_num = list(map(lambda x: x - 1, mask_num))  # Convertiamo in indici 0-based
                mask_out = list(map(lambda x: x + 1, mask_num))  # Per la visualizzazione
            except:
                print("Errore nella conversione della maschera personalizzata.")
                mask_num = [0, 1, -1]  # Default: prima, seconda e ultima posizione
                mask_out = [1, 2, "cterm"]
        
        # Stampiamo le informazioni sulla predizione
        if allele:
            print(f"Allele: {allele}")
        print(f"Mascheramento: {'personalizzato' if custom_mask else 'predefinito'}")
        print(f"Posizioni mascherate: {mask_out}\n")
        
        # Processiamo ogni peptide
        for peptide in tqdm(valid_peptides, desc="Predizione immunogenicità"):
            peptide = peptide.upper()
            peplen = len(peptide)
            
            cterm = peplen - 1
            score = 0
            count = 0
            
            # Determiniamo la maschera effettiva
            if not custom_mask:
                mask_num = [0, 1, cterm]
            
            # Adattiamo i pesi in base alla lunghezza del peptide
            if peplen > 9:
                pepweight = self.immunoweight[:5] + ((peplen - 9) * [0.30]) + self.immunoweight[5:]
            else:
                pepweight = self.immunoweight
            
            try:
                # Calcoliamo lo score
                for pos, aa in enumerate(peptide):
                    if aa not in self.immunoscale:
                        raise KeyError(f"Amminoacido non valido: {aa}")
                    elif pos not in mask_num:
                        score += pepweight[pos] * self.immunoscale[aa]
                
                # Aggiungiamo il risultato alla lista
                result_list.append({
                    'peptide': peptide,
                    'length': peplen,
                    'score': round(score, 5)
                })
            
            except Exception as e:
                print(f"Errore nell'analisi del peptide {peptide}: {e}")
        
        # Ordiniamo i risultati per score decrescente
        result_list.sort(key=lambda x: x['score'], reverse=True)
        
        return result_list
    
    def save_results_to_csv(self, results, filename="immunogenicity_results.csv"):
        """
        Salva i risultati della predizione in un file CSV.
        
        Parametri:
        results (list): Lista di risultati della predizione
        filename (str): Nome del file CSV
        
        Returns:
        str: Nome del file CSV creato
        """
        if not results:
            print("Nessun risultato da salvare.")
            return None
        
        try:
            with open(filename, 'w', newline='') as csvfile:
                fieldnames = ['peptide', 'length', 'score']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for result in results:
                    writer.writerow(result)
            
            print(f"Risultati salvati in {filename}")
            return filename
        
        except Exception as e:
            print(f"Errore nel salvataggio dei risultati in {filename}: {e}")
            return None
    
    def get_available_alleles(self):
        """
        Restituisce la lista degli alleli disponibili.
        
        Returns:
        list: Lista degli alleli disponibili
        """
        return sorted(list(self.allele_dict.keys()))
    
    def print_available_alleles(self):
        """
        Stampa la lista degli alleli disponibili in formato tabellare.
        """
        alleles = self.get_available_alleles()
        
        print("\n========================================================================= ")
        print("| Alleli disponibili per la predizione dell'immunogenicità di Classe I:  | ")
        print("|-----------------------------------------------------------------------| ")
        
        # Stampiamo gli alleli in righe da 6
        for i in range(0, len(alleles), 6):
            row = alleles[i:i+6]
            print("| " + " | ".join(f"{allele:9}" for allele in row) + " | ")
        
        print("-------------------------------------------------------------------------\n")

def predict_peptide_immunogenicity(peptides, custom_mask=None, allele=None, output_csv=None):
    """
    Funzione di utilità per predire l'immunogenicità di una lista di peptidi.
    
    Parametri:
    peptides (list): Lista di peptidi da analizzare
    custom_mask (str): Maschera personalizzata (opzionale)
    allele (str): Allele HLA da utilizzare (opzionale)
    output_csv (str): Nome del file CSV di output (opzionale)
    
    Returns:
    list: Lista di risultati della predizione
    """
    predictor = Immunogenicity()
    results = predictor.predict_immunogenicity(peptides, custom_mask, allele)
    
    if output_csv and results:
        predictor.save_results_to_csv(results, output_csv)
    
    return results

def get_available_immunogenicity_alleles():
    """
    Restituisce la lista degli alleli disponibili per la predizione dell'immunogenicità.
    
    Returns:
    list: Lista degli alleli disponibili
    """
    predictor = Immunogenicity()
    return predictor.get_available_alleles()

def print_available_immunogenicity_alleles():
    """
    Stampa la lista degli alleli disponibili per la predizione dell'immunogenicità.
    """
    predictor = Immunogenicity()
    predictor.print_available_alleles()
