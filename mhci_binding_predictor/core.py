import os
import logging
import tempfile
import numpy as np
import pandas as pd
import requests
import time

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("IEDBBindingPredictor")

class IEDBBindingPredictor:
    """
    Class for predicting peptide binding with MHC alleles using the IEDB API.
    Also calculates immunogenicity scores based on amino acid properties.
    """
    
    def __init__(self, output_dir="./output", method="recommended_epitope", csv_separator=",", decimal_separator="."):
        """
        Initialize the binding predictor.
        
        Args:
            output_dir (str): Output directory for results
            method (str): Prediction method to use (can include version with dash, e.g. "netmhcpan_el-4.1")
            csv_separator (str): Separator to use for CSV files ("," or ";")
            decimal_separator (str): Decimal separator to use for CSV files ("." or ",")
        """
        self.output_dir = output_dir
        
        # Gestione versione metodo
        if "-" in method:
            self.method, self.method_version = method.split("-", 1)
        else:
            self.method = method
            self.method_version = None
        
        # Set CSV export settings
        if csv_separator in [",", ";"]:
            self.csv_separator = csv_separator
        else:
            logger.warning(f"Invalid CSV separator '{csv_separator}', using ',' instead")
            self.csv_separator = ","
            
        if decimal_separator in [".", ","]:
            self.decimal_separator = decimal_separator
        else:
            logger.warning(f"Invalid decimal separator '{decimal_separator}', using '.' instead")
            self.decimal_separator = "."
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize immunogenicity calculation variables
        self.immunoscale = {
            "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110, 
            "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036, "M": -0.570, "N": -0.021, 
            "P": -0.036, "Q": -0.376, "R": 0.168, "S": -0.537, "T": 0.126, "V": 0.134, 
            "W": 0.719, "Y": -0.012
        }
        self.immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
        
        # Initialize allele dictionary for masks
        self.allele_dict = {
            "H-2-Db":"2,5,9", "H-2-Dd":"2,3,5", "H-2-Kb":"2,3,9", "H-2-Kd":"2,5,9", 
            "H-2-Kk":"2,8,9", "H-2-Ld":"2,5,9", "HLA-A0101":"2,3,9", "HLA-A0201":"1,2,9", 
            "HLA-A0202":"1,2,9", "HLA-A0203":"1,2,9", "HLA-A0206":"1,2,9", "HLA-A0211":"1,2,9", 
            "HLA-A0301":"1,2,9", "HLA-A1101":"1,2,9", "HLA-A2301":"2,7,9", "HLA-A2402":"2,7,9", 
            "HLA-A2601":"1,2,9", "HLA-A2902":"2,7,9", "HLA-A3001":"1,3,9", "HLA-A3002":"2,7,9", 
            "HLA-A3101":"1,2,9", "HLA-A3201":"1,2,9", "HLA-A3301":"1,2,9", "HLA-A6801":"1,2,9", 
            "HLA-A6802":"1,2,9", "HLA-A6901":"1,2,9", "HLA-B0702":"1,2,9", "HLA-B0801":"2,5,9", 
            "HLA-B1501":"1,2,9", "HLA-B1502":"1,2,9", "HLA-B1801":"1,2,9", "HLA-B2705":"2,3,9", 
            "HLA-B3501":"1,2,9", "HLA-B3901":"1,2,9", "HLA-B4001":"1,2,9", "HLA-B4002":"1,2,9", 
            "HLA-B4402":"2,3,9", "HLA-B4403":"2,3,9", "HLA-B4501":"1,2,9", "HLA-B4601":"1,2,9", 
            "HLA-B5101":"1,2,9", "HLA-B5301":"1,2,9", "HLA-B5401":"1,2,9", "HLA-B5701":"1,2,9",
            "HLA-B5801":"1,2,9"
        }
        
        logger.info(f"Initialized binding predictor with method: {method}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"CSV separator: '{self.csv_separator}', Decimal separator: '{self.decimal_separator}'")
    
    def isint(self, value):
        """
        Check if a value is an integer.
        
        Args:
            value: Value to check
            
        Returns:
            bool: True if the value is an integer, False otherwise
        """
        try:
            int(value)
            return True
        except ValueError:
            return False
    
    def predict_binding_api(self, peptides, alleles, lengths=None, method=None):
        """
        Predicts peptide binding with the specified alleles using the IEDB API.
        Supporta richieste batch con multipli alleli e lunghezze.
        
        Args:
            peptides (list): List of peptides
            alleles (list or str): MHC allele(s)
            lengths (list or int, optional): Peptide length(s)
            method (str, optional): Override the default method
            
        Returns:
            pd.DataFrame: DataFrame with prediction results
        """

        time.sleep(1)

        if method is None:
            method = self.method
            
        # Converti parametri a liste se necessario
        if isinstance(alleles, str):
            alleles = [alleles]
            
        if lengths is None:
            lengths = [9]
        elif isinstance(lengths, int):
            lengths = [lengths]
            
        logger.info(f"Using IEDB API for binding prediction with {len(alleles)} alleles and method {method}")
        
        # IEDB API URL
        url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
        
        # Costruisci parametri versione
        method_param = method
        if self.method_version and method in ["netmhcpan_el", "netmhcpan_ba"]:
            method_param = f"{method}-{self.method_version}"
            
        # Formatta correttamente sequenze multiple
        fasta_sequences = "\n".join([f">peptide{i}\n{p}" for i, p in enumerate(peptides)])
        
        # Prepare data for the request
        data = {
            "method": method_param,
            "sequence_text": fasta_sequences,
            "allele": ",".join(alleles),
            "length": ",".join(map(str, lengths))
        }
        
        try:
            # Make POST request
            response = requests.post(url, data=data)
            
            # Verify request was successful
            if response.status_code == 200:
                # Convert response to DataFrame
                results = response.text.strip().split("\n")
                
                # Verify there are results
                if len(results) < 2:
                    logger.warning(f"No results from IEDB API for allele {allele}")
                    return pd.DataFrame()
                
                # Create temporary file with results
                with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.csv') as temp_file:
                    temp_file.write(response.text)
                    temp_file_path = temp_file.name
                
                # Read CSV file
                df = pd.read_csv(temp_file_path, sep="\t")
                
                # Remove temporary file
                os.remove(temp_file_path)
                
                # Mantieni nomi colonna originali e seleziona solo colonne rilevanti
                relevant_columns = ['allele', 'peptide', 'start', 'end', 'method', 'percentile_rank']
                
                # Aggiungi colonne specifiche del metodo se presenti
                for col in ['ann_ic50', 'netmhcpan_ic50', 'netmhcpan_rank', 
                           'netmhcpan_el_score', 'netmhcpan_el_rank', 
                           'netmhcpan_ba_score', 'netmhcpan_ba_rank', 'netmhcpan_ba_ic50']:
                    if col in df.columns:
                        relevant_columns.append(col)
                
                # Seleziona solo le colonne rilevanti se presenti
                df = df[[col for col in relevant_columns if col in df.columns]]
                
                # Save results to CSV file
                api_csv = os.path.join(self.output_dir, f"api_binding_{allele.replace('*', '').replace(':', '')}_{method}.csv")
                self.save_to_csv(df, api_csv)
                logger.info(f"API results saved to: {api_csv}")
                
                return df
            else:
                logger.error(f"Error in IEDB API request: {response.status_code} - {response.text}")
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"Error parsing API response: {str(e)}")
            # Tentare un approccio alternativo per il parsing
            try:
                lines = response.text.strip().split("\n")
                if len(lines) >= 2:
                    headers = lines[0].split("\t")
                    data = []
                    for line in lines[1:]:
                        values = line.split("\t")
                        if len(values) == len(headers):
                            data.append(dict(zip(headers, values)))
                    df = pd.DataFrame(data)
                    return df
            except:
                return pd.DataFrame()
    
    def calculate_immunogenicity_score(self, peptide, custom_mask=None, allele=None, method=None):
        """
        Calculate immunogenicity score for a peptide.
        
        Args:
            peptide (str): Peptide sequence
            custom_mask (str, optional): Custom mask positions (comma-separated)
            allele (str, optional): MHC allele name
            method (str, optional): Prediction method used
            
        Returns:
            float: Immunogenicity score
        """
        # Logica adattiva basata sul metodo
        if method and "netmhcpan_el" in method:
            return self._calculate_modern_immunogenicity(peptide, allele, custom_mask)
        else:
            return self._calculate_legacy_immunogenicity(peptide, allele, custom_mask)
    
    def _calculate_legacy_immunogenicity(self, peptide, allele=None, custom_mask=None):
        """
        Calcolo immunogenicità con algoritmo legacy.
        """
        peptide = peptide.upper()
        peplen = len(peptide)
        cterm = peplen - 1
        score = 0
        count = 0
        
        # Set mask based on allele or custom mask or default
        if allele and allele in self.allele_dict:
            mask_str = self.allele_dict[allele].split(",")
            mask_num = list(map(int, mask_str))
            mask_num = list(map(lambda x: x - 1, mask_num))  # Convert to 0-based indexing
        elif custom_mask:
            mask_str = custom_mask.split(",")
            mask_num = list(map(int, mask_str))
            mask_num = list(map(lambda x: x - 1, mask_num))  # Convert to 0-based indexing
        else:
            mask_num = [0, 1, cterm]  # Default mask (first, second, and last positions)
        
        # Adjust weights for peptides longer than 9 amino acids
        if peplen > 9:
            pepweight = self.immunoweight[:5] + ((peplen - 9) * [0.30]) + self.immunoweight[5:]
        else:
            pepweight = self.immunoweight
        
        try:
            for count, aa in enumerate(peptide):
                if aa not in self.immunoscale:
                    logger.warning(f"Invalid amino acid '{aa}' in peptide {peptide}")
                    return 0.0  # Return 0 score for invalid amino acids
                elif count not in mask_num:
                    score += pepweight[min(count, len(pepweight)-1)] * self.immunoscale[aa]
            
            return round(score, 5)
        
        except Exception as e:
            logger.error(f"Error calculating immunogenicity for {peptide}: {str(e)}")
            return 0.0
    
    def _calculate_modern_immunogenicity(self, peptide, allele=None, custom_mask=None):
        """
        Calcolo immunogenicità con algoritmo moderno ottimizzato per netmhcpan_el.
        Utilizza pesi e scale aggiornate per una migliore correlazione con i dati sperimentali.
        
        Args:
            peptide (str): Peptide sequence
            allele (str, optional): MHC allele name
            custom_mask (str, optional): Custom mask positions
            
        Returns:
            float: Immunogenicity score
        """
        peptide = peptide.upper()
        peplen = len(peptide)
        
        # Scala immunogenicità moderna (2023)
        modern_scale = {
            "A": 0.107, "C": -0.205, "D": 0.102, "E": 0.355, "F": 0.410, 
            "G": 0.130, "H": 0.125, "I": 0.462, "K": -0.730, "L": -0.016, 
            "M": -0.550, "N": -0.001, "P": -0.056, "Q": -0.356, "R": 0.188, 
            "S": -0.517, "T": 0.146, "V": 0.154, "W": 0.739, "Y": 0.008
        }
        
        # Pesi posizionali moderni con maggiore enfasi sulle posizioni centrali
        modern_weights = [0.05, 0.10, 0.25, 0.40, 0.45, 0.45, 0.40, 0.25, 0.10]
        
        # Estendi i pesi per peptidi più lunghi di 9 aa
        if peplen > 9:
            middle_weights = [0.45] * (peplen - 9)
            weights = modern_weights[:4] + middle_weights + modern_weights[5:]
        else:
            weights = modern_weights[:peplen]
        
        # Determina le posizioni di ancoraggio da mascherare
        if allele and allele.replace("*", "").replace(":", "") in self.allele_dict:
            formatted_allele = allele.replace("*", "").replace(":", "")
            mask_str = self.allele_dict[formatted_allele].split(",")
            mask_positions = list(map(int, mask_str))
            mask_positions = list(map(lambda x: x - 1, mask_positions))  # 0-based
        elif custom_mask:
            mask_str = custom_mask.split(",")
            mask_positions = list(map(int, mask_str))
            mask_positions = list(map(lambda x: x - 1, mask_positions))  # 0-based
        else:
            # Maschera predefinita migliorata
            mask_positions = [0, 1, peplen-1]
        
        try:
            score = 0.0
            for i, aa in enumerate(peptide):
                if aa not in modern_scale:
                    logger.warning(f"Amminoacido non valido '{aa}' nel peptide {peptide}")
                    continue
                
                # Applica il punteggio solo per posizioni non di ancoraggio
                if i not in mask_positions:
                    position_weight = weights[min(i, len(weights)-1)]
                    score += position_weight * modern_scale[aa]
            
            # Normalizza il punteggio in base alla lunghezza
            normalized_score = score / (peplen - len(mask_positions))
            return round(normalized_score * 2, 5)  # Fattore di scala per allineare con i punteggi legacy
            
        except Exception as e:
            logger.error(f"Errore nel calcolo dell'immunogenicità moderna per {peptide}: {str(e)}")
            return 0.0
    
    def add_immunogenicity(self, df):
        """
        Add immunogenicity score column to DataFrame.
        """
        if df.empty or 'peptide' not in df.columns:
            return df
        
        try:
            # Calcola l'immunogenicità per ogni peptide usando allele e metodo specifici
            df['immunogenicity'] = df.apply(
                lambda row: self.calculate_immunogenicity_score(
                    row['peptide'],
                    allele=row['allele'].replace("*", "").replace(":", "") if 'allele' in row else None,
                    method=row['method'] if 'method' in row else None
                ), 
                axis=1
            )
            
            logger.info(f"Added immunogenicity scores for {len(df)} peptides")
            return df
        
        except Exception as e:
            logger.error(f"Error adding immunogenicity: {str(e)}")
            df['immunogenicity'] = np.nan  # Aggiungi colonna vuota in caso di errore
            return df
    
    def save_to_csv(self, df, file_path):
        """
        Save DataFrame to CSV with custom separator and decimal separator.
        
        Args:
            df (pd.DataFrame): DataFrame to save
            file_path (str): Path where to save the CSV file
        """
        if df.empty:
            logger.warning(f"Empty DataFrame, not saving to {file_path}")
            return
            
        try:
            # Using the configured separators
            df.to_csv(file_path, 
                     sep=self.csv_separator, 
                     decimal=self.decimal_separator, 
                     index=False)
            logger.info(f"Saved DataFrame to {file_path} with separator '{self.csv_separator}' and decimal '{self.decimal_separator}'")
        except Exception as e:
            logger.error(f"Error saving to CSV: {str(e)}")
    
    def analyze_peptides(self, peptides, alleles, lengths=None):
        """
        Analyze peptides using the configured method and return results.
        Supporta richieste batch con multipli alleli e lunghezze.
        
        Args:
            peptides (list): List of peptides to analyze
            alleles (list or str): MHC allele(s)
            lengths (list or int, optional): Peptide length(s)
            
        Returns:
            pd.DataFrame: DataFrame with analysis results
        """
        # Converti parametri a liste se necessario
        if isinstance(alleles, str):
            alleles = [alleles]
            
        if lengths is None:
            lengths = [9]
        elif isinstance(lengths, int):
            lengths = [lengths]
            
        logger.info(f"Analyzing {len(peptides)} peptides with {len(alleles)} alleles using method {self.method}")
        
        # Ottieni predizioni in un'unica chiamata API
        results = self.predict_binding_api(peptides, alleles, lengths)
        
        if results.empty:
            logger.error("Failed to get predictions")
            return pd.DataFrame()
        
        try:
            # Aggiungi immunogenicità
            results = self.add_immunogenicity(results)
            
            # Mappatura colonne dinamica in base al metodo
            if 'netmhcpan_el_score' in results.columns and 'el' in self.method:
                results['score'] = results['netmhcpan_el_score']
            elif 'netmhcpan_ba_ic50' in results.columns and 'ba' in self.method:
                results['score'] = results['netmhcpan_ba_ic50']
            elif 'ann_ic50' in results.columns:
                results['score'] = results['ann_ic50']
            
            # Aggiungi colonna IC50 se disponibile
            if 'netmhcpan_ba_ic50' in results.columns:
                results['ic50'] = results['netmhcpan_ba_ic50']
            elif 'netmhcpan_ic50' in results.columns:
                results['ic50'] = results['netmhcpan_ic50']
            elif 'ann_ic50' in results.columns:
                results['ic50'] = results['ann_ic50']
            else:
                results['ic50'] = np.nan
            
            # Assicurati che tutte le colonne richieste siano presenti
            if 'percentile_rank' not in results.columns:
                results['percentile_rank'] = np.nan
            if 'score' not in results.columns:
                results['score'] = np.nan
                
            # Seleziona e ordina le colonne finali
            final_columns = ['peptide', 'allele', 'score', 'percentile_rank', 'immunogenicity', 'ic50']
            final_df = results[final_columns]
            
            # Salva i risultati
            results_csv = os.path.join(self.output_dir, f"analysis_{self.method}.csv")
            self.save_to_csv(final_df, results_csv)
            logger.info(f"Analysis results saved to: {results_csv}")
            
            return final_df
            
        except Exception as e:
            logger.error(f"Error merging predictions: {str(e)}")
            # Return one of the results sets if available
            if not el_results.empty:
                logger.info("Returning only EL results due to merge error")
                return el_results
            elif not ba_results.empty:
                logger.info("Returning only BA results due to merge error")
                return ba_results
            else:
                return pd.DataFrame()
    
    def get_required_columns(self, df):
        """
        Ensure consistent output columns.
        
        Args:
            df (pd.DataFrame): DataFrame to standardize
            
        Returns:
            pd.DataFrame: DataFrame with required columns
        """
        required_cols = ['peptide', 'allele', 'score', 'percentile_rank', 'immunogenicity', 'ic50']
        for col in required_cols:
            if col not in df.columns:
                df[col] = np.nan
        return df[required_cols]
    
    def analyze_peptides_batch(self, peptides_list, alleles, lengths=None, batch_size=50):
        """
        Analyze batches of peptides across multiple alleles.
        
        Args:
            peptides_list (list): List of peptides
            alleles (list or str): MHC allele(s)
            lengths (list or int, optional): Peptide length(s)
            batch_size (int): Batch size
            
        Returns:
            pd.DataFrame: Combined DataFrame with results for all alleles
        """
        # Converti parametri a liste se necessario
        if isinstance(alleles, str):
            alleles = [alleles]
            
        if lengths is None:
            lengths = [9]
        elif isinstance(lengths, int):
            lengths = [lengths]
        
        all_results = []
        
        # Split peptides into batches to avoid overloading the API
        for i in range(0, len(peptides_list), batch_size):
            batch_peptides = peptides_list[i:i+batch_size]
            logger.info(f"Processing batch {i//batch_size + 1} with {len(batch_peptides)} peptides")
            
            # Analizza tutti gli alleli e lunghezze in un'unica chiamata per batch
            batch_result = self.analyze_peptides(batch_peptides, alleles, lengths)
            if not batch_result.empty:
                all_results.append(batch_result)
        
        # Combina tutti i risultati in un unico DataFrame
        if all_results:
            combined_results = pd.concat(all_results, ignore_index=True)
            combined_results = self.get_required_columns(combined_results)
            
            # Salva i risultati combinati
            combined_csv = os.path.join(self.output_dir, f"all_batches_{self.method}.csv")
            self.save_to_csv(combined_results, combined_csv)
            logger.info(f"Combined {len(combined_results)} results from all batches")
            
            return combined_results
        else:
            logger.warning("No results to combine")
            return pd.DataFrame(columns=['peptide', 'allele', 'score', 'percentile_rank', 'immunogenicity', 'ic50'])
    
    def filter_binders(self, df, threshold=0., rank_threshold=None, ic50_threshold=None, immunogenicity_threshold=None):
        """
        Filter binders based on score, rank, IC50 or immunogenicity.
        
        Args:
            df (pd.DataFrame): DataFrame with prediction results
            threshold (float): Threshold for score (higher = better binding)
            rank_threshold (float): Threshold for percentile rank (lower = better binding)
            ic50_threshold (float): Threshold for IC50 (lower = better binding)
            immunogenicity_threshold (float, optional): Threshold for immunogenicity (higher = more immunogenic)
            
        Returns:
            pd.DataFrame: DataFrame with filtered binders
        """
        if df.empty:
            return df

        # Converti le colonne numeriche usando il separatore decimale corrente
        numeric_cols = ['score', 'percentile_rank', 'ic50', 'immunogenicity']
        for col in numeric_cols:
            if col in df.columns and df[col].dtype == object:
                # Sostituisci il separatore decimale con un punto e converti a float
                df[col] = df[col].str.replace(self.decimal_separator, '.', regex=False)
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Filter based on available metrics
        filters = []
        
        if 'score' in df.columns and threshold is not None:
            filters.append(df['score'] >= threshold)
        
        if 'percentile_rank' in df.columns and rank_threshold is not None:
            filters.append(df['percentile_rank'] <= rank_threshold)
        
        if 'ic50' in df.columns and ic50_threshold is not None:
            filters.append(df['ic50'] <= ic50_threshold)
            
        if 'immunogenicity' in df.columns and immunogenicity_threshold is not None:
            filters.append(df['immunogenicity'] >= immunogenicity_threshold)
            
        # Apply filters if any exist
        if filters:
            # Combine all filters with AND (changed from OR to AND for more restrictive filtering)
            combined_filter = filters[0]
            for f in filters[1:]:
                combined_filter &= f
            
            binders = df[combined_filter].copy()
            logger.info(f"Identificati {len(binders)} binder che soddisfano tutte le soglie")
            return binders
        else:
            logger.warning("Impossibile filtrare: nessuna colonna valida per il filtraggio")
            return df
    
    def set_csv_format(self, csv_separator=None, decimal_separator=None):
        """
        Update CSV formatting options.
        
        Args:
            csv_separator (str, optional): Separator for CSV files ("," or ";")
            decimal_separator (str, optional): Decimal separator for CSV files ("." or ",")
        """
        if csv_separator is not None:
            if csv_separator in [",", ";"]:
                self.csv_separator = csv_separator
                logger.info(f"CSV separator set to '{csv_separator}'")
            else:
                logger.warning(f"Invalid CSV separator '{csv_separator}', keeping current setting")
        
        if decimal_separator is not None:
            if decimal_separator in [".", ","]:
                self.decimal_separator = decimal_separator
                logger.info(f"Decimal separator set to '{decimal_separator}'")
            else:
                logger.warning(f"Invalid decimal separator '{decimal_separator}', keeping current setting")

    def generate_variants(self, sequence_pattern):
        """
        Generate variants from a pattern (e.g. A[CD]E[FY]GH).
        
        Args:
            sequence_pattern (str): Pattern with options in brackets
            
        Returns:
            list: List of all possible peptide variants
        """
        variants = []
        
        # Recursive function to generate all combinations
        def generate_combinations(current_seq, remaining_pattern):
            if not remaining_pattern:
                variants.append(current_seq)
                return
                
            # If character starts with '[', extract options between brackets
            if remaining_pattern[0] == '[':
                close_bracket = remaining_pattern.find(']')
                if close_bracket > 0:
                    options = remaining_pattern[1:close_bracket]
                    for option in options:
                        generate_combinations(
                            current_seq + option,
                            remaining_pattern[close_bracket+1:]
                        )
                    return
            
            # Normal character
            generate_combinations(
                current_seq + remaining_pattern[0],
                remaining_pattern[1:]
            )
                
        generate_combinations("", sequence_pattern)
        
        if not variants:
            # If no variants, consider the pattern as direct sequence
            variants = [sequence_pattern]
            
        logger.info(f"Generated {len(variants)} variants from pattern sequence")
        return variants

def main():
    """Main entry point for CLI"""
    pass
