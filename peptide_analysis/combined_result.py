"""
Modulo per la gestione dei risultati combinati di binding MHC-I e immunogenicità.
"""

import csv
import pandas as pd
from .immunogenicity import predict_peptide_immunogenicity

class CombinedResult:
    """
    Classe per la gestione dei risultati combinati di binding MHC-I e immunogenicità.
    
    Questa classe fornisce metodi per calcolare il punteggio composito e categorizzare
    i peptidi in base ai loro valori di binding MHC-I e immunogenicità.
    """
    
    def __init__(self, peptide, allele, binding_score, percentile_rank, immunogenicity_score=None):
        """
        Inizializza un oggetto CombinedResult.
        
        Parametri:
        peptide (str): Sequenza peptidica
        allele (str): Allele HLA
        binding_score (float): Score di binding MHC-I
        percentile_rank (float): Percentile rank del binding MHC-I
        immunogenicity_score (float): Score di immunogenicità (opzionale)
        """
        self.peptide = peptide
        self.allele = allele
        self.binding_score = binding_score
        self.percentile_rank = percentile_rank
        self.immunogenicity_score = immunogenicity_score
        self.composite_score = self._calculate_composite_score()
        self.category = self._determine_category()
    
    def _calculate_composite_score(self):
        """
        Calcola il punteggio composito che combina binding MHC-I e immunogenicità.
        
        La formula utilizzata è:
        Punteggio composito = (immunogenicity_score * 0.5) + ((1 - percentile_rank/100) * 0.3) + (binding_score * 0.2)
        
        Returns:
        float: Punteggio composito
        """
        if self.immunogenicity_score is None:
            return None
        
        # Formula del punteggio composito
        immunogenicity_component = self.immunogenicity_score * 0.5
        binding_rank_component = (1 - (self.percentile_rank / 100)) * 0.3
        binding_score_component = self.binding_score * 0.2
        
        # Il punteggio è una media ponderata dei tre componenti
        composite_score = immunogenicity_component + binding_rank_component + binding_score_component
        
        return round(composite_score, 4)
    
    def _determine_category(self):
        """
        Determina la categoria del peptide in base al percentile rank e all'immunogenicità.
        
        Returns:
        str: Categoria del peptide
        """
        if self.immunogenicity_score is None:
            return 'Non determinato'
        
        # Categorizzazione in base ai criteri definiti
        if self.immunogenicity_score > 0.3 and self.percentile_rank < 0.1 and self.binding_score > 0.95:
            return 'Eccellente'
        elif self.immunogenicity_score > 0 and self.percentile_rank < 0.5 and self.binding_score > 0.9:
            return 'Buono'
        elif self.immunogenicity_score > 0 and self.percentile_rank < 1.0 and self.binding_score > 0.8:
            return 'Da considerare'
        else:
            return 'Da scartare'
    
    def to_dict(self):
        """
        Converte l'oggetto in un dizionario.
        
        Returns:
        dict: Dizionario con i dati del risultato
        """
        return {
            'peptide': self.peptide,
            'allele': self.allele,
            'score': self.binding_score,
            'percentile_rank': self.percentile_rank,
            'immunogenicity_score': self.immunogenicity_score,
            'punteggio_composito': self.composite_score,
            'categoria': self.category
        }
    
    @classmethod
    def from_dict(cls, data):
        """
        Crea un oggetto CombinedResult da un dizionario.
        
        Parametri:
        data (dict): Dizionario con i dati del risultato
        
        Returns:
        CombinedResult: Oggetto CombinedResult
        """
        return cls(
            peptide=data['peptide'],
            allele=data['allele'],
            binding_score=data['score'],
            percentile_rank=data['percentile_rank'],
            immunogenicity_score=data.get('immunogenicity_score')
        )
    
    @classmethod
    def from_results_list(cls, results_list):
        """
        Crea una lista di oggetti CombinedResult da una lista di dizionari.
        
        Parametri:
        results_list (list): Lista di dizionari con i dati dei risultati
        
        Returns:
        list: Lista di oggetti CombinedResult
        """
        return [cls.from_dict(result) for result in results_list]
    
    @staticmethod
    def add_immunogenicity_scores(results, output_csv=None, filtered_csv=None, ranked_csv=None):
        """
        Aggiunge il punteggio di immunogenicità ai risultati dell'analisi.
        
        Parametri:
        results (list): Lista di dizionari con i risultati dell'analisi
        output_csv (str): Nome del file CSV di output (opzionale)
        filtered_csv (str): Nome del file CSV filtrato (opzionale)
        ranked_csv (str): Nome del file CSV con i risultati ordinati (opzionale)
        
        Returns:
        tuple: (risultati con immunogenicità, risultati filtrati, risultati ordinati)
        """
        # Estrai i peptidi unici dai risultati
        unique_peptides = list(set([r['peptide'] for r in results]))
        
        # Analizziamo l'immunogenicità dei peptidi
        immuno_results = predict_peptide_immunogenicity(unique_peptides)
        
        # Creiamo un dizionario per mappare peptidi a score di immunogenicità
        immuno_scores = {r['peptide']: r['score'] for r in immuno_results}
        
        # Aggiungiamo lo score di immunogenicità ai risultati originali
        combined_results = []
        for r in results:
            peptide = r['peptide']
            allele = r['allele']
            allele_clean = allele.replace('*', '').replace(':', '')
            
            # Calcoliamo l'immunogenicità specifica per l'allele
            allele_specific_results = predict_peptide_immunogenicity([peptide], allele=allele_clean)
            if allele_specific_results:
                immunogenicity_score = allele_specific_results[0]['score']
            else:
                # Fallback al valore generico se non disponibile per l'allele specifico
                immunogenicity_score = immuno_scores.get(peptide, None)
            
            # Creiamo un oggetto CombinedResult
            combined_result = CombinedResult(
                peptide=peptide,
                allele=allele,
                binding_score=r['score'],
                percentile_rank=r['percentile_rank'],
                immunogenicity_score=immunogenicity_score
            )
            
            combined_results.append(combined_result)
        
        # Convertiamo gli oggetti CombinedResult in dizionari
        results_with_immuno = [cr.to_dict() for cr in combined_results]
        
        # Filtriamo i risultati (percentile_rank < 0.5 e immunogenicity > 0)
        filtered_results = [
            cr.to_dict() for cr in combined_results 
            if cr.immunogenicity_score is not None and cr.immunogenicity_score > 0 and cr.percentile_rank < 0.5
        ]
        
        # Ordiniamo i risultati per punteggio composito decrescente
        ranked_results = sorted(
            [cr.to_dict() for cr in combined_results if cr.composite_score is not None],
            key=lambda x: x['punteggio_composito'],
            reverse=True
        )
        
        # Salviamo i risultati in file CSV se richiesto
        if output_csv and results_with_immuno:
            df = pd.DataFrame(results_with_immuno)
            df.to_csv(output_csv, index=False)
            print(f"File CSV con punteggi di immunogenicità creato: {output_csv}")
        
        if filtered_csv and filtered_results:
            df = pd.DataFrame(filtered_results)
            df.to_csv(filtered_csv, index=False)
            print(f"File CSV filtrato creato: {filtered_csv}")
            print(f"Numero di peptidi nel file filtrato: {len(filtered_results)}")
        
        if ranked_csv and ranked_results:
            df = pd.DataFrame(ranked_results)
            df.to_csv(ranked_csv, index=False)
            print(f"File CSV con peptidi classificati creato: {ranked_csv}")
        
        return results_with_immuno, filtered_results, ranked_results
