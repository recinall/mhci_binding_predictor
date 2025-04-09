import os
import logging
import tempfile
import numpy as np
import pandas as pd
import requests
import sys

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
    
    def __init__(self, output_dir="./output", method="netmhcpan_el", csv_separator=",", decimal_separator="."):
        """
        Initialize the binding predictor.
        
        Args:
            output_dir (str): Output directory for results
            method (str): Prediction method to use
            csv_separator (str): Separator to use for CSV files ("," or ";")
            decimal_separator (str): Decimal separator to use for CSV files ("." or ",")
        """
        self.output_dir = output_dir
        self.method = method
        
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
    
    def predict_binding_api(self, peptides, allele, length=9, method=None):
        """
        Predicts peptide binding with the specified allele using the IEDB API.
        
        Args:
            peptides (list): List of peptides
            allele (str): MHC allele
            length (int): Peptide length
            method (str, optional): Override the default method
            
        Returns:
            pd.DataFrame: DataFrame with prediction results
        """
        if method is None:
            method = self.method
            
        logger.info(f"Using IEDB API for binding prediction with allele {allele} and method {method}")
        
        # IEDB API URL
        url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
        
        # Prepare data for the request
        data = {
            "method": method if method else "recommended",
            "sequence_text": "\n".join(peptides),
            "allele": allele,
            "length": str(length)
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
                
                # Ensure all required metrics are present
                # Rename columns if necessary to standardize names
                column_mapping = {
                    'start': 'start',
                    'end': 'end',
                    'peptide': 'peptide',
                    'method': 'method',
                    'percentile_rank': 'percentile_rank',
                    'ann_ic50': 'ann_ic50',
                    'ann_rank': 'ann_rank',
                    'smm_ic50': 'smm_ic50',
                    'smm_rank': 'smm_rank',
                    'comblib_sidney2008_score': 'comblib_score',
                    'comblib_sidney2008_rank': 'comblib_rank',
                    'netmhcpan_ic50': 'netmhcpan_ic50',
                    'netmhcpan_rank': 'netmhcpan_rank',
                    'netmhcpan_el_score': 'netmhcpan_el_score',
                    'netmhcpan_el_rank': 'netmhcpan_el_rank',
                    'netmhcpan_ba_score': 'netmhcpan_ba_score',
                    'netmhcpan_ba_rank': 'netmhcpan_ba_rank',
                    'netmhcpan_ba_ic50': 'netmhcpan_ba_ic50'
                }
                
                # Rename columns if present
                for old_col, new_col in column_mapping.items():
                    if old_col in df.columns and new_col not in df.columns:
                        df[new_col] = df[old_col]
                
                # Add allele information
                df['allele'] = allele
                
                # Save results to CSV file
                api_csv = os.path.join(self.output_dir, f"api_binding_{allele.replace('*', '').replace(':', '')}_{method}.csv")
                self.save_to_csv(df, api_csv)
                logger.info(f"API results saved to: {api_csv}")
                
                return df
            else:
                logger.error(f"Error in IEDB API request: {response.status_code} - {response.text}")
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"Error during prediction with IEDB API: {str(e)}")
            return pd.DataFrame()
    
    def calculate_immunogenicity_score(self, peptide, custom_mask=None, allele=None):
        """
        Calculate immunogenicity score for a peptide.
        
        Args:
            peptide (str): Peptide sequence
            custom_mask (str, optional): Custom mask positions (comma-separated)
            allele (str, optional): MHC allele name
            
        Returns:
            float: Immunogenicity score
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
    
    def add_immunogenicity(self, df, allele=None):
        """
        Add immunogenicity score column to DataFrame.
        
        Args:
            df (pd.DataFrame): DataFrame with peptides
            allele (str, optional): MHC allele name
            
        Returns:
            pd.DataFrame: DataFrame with added immunogenicity scores
        """
        if df.empty or 'peptide' not in df.columns:
            return df
        
        try:
            # Calculate immunogenicity for each peptide
            df['immunogenicity'] = df['peptide'].apply(
                lambda x: self.calculate_immunogenicity_score(x, allele=allele)
            )
            
            logger.info(f"Added immunogenicity scores for {len(df)} peptides")
            return df
        
        except Exception as e:
            logger.error(f"Error adding immunogenicity: {str(e)}")
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
    
    def analyze_peptides(self, peptides, allele, length=9):
        """
        Analyze peptides using both netmhcpan_el and netmhcpan_ba methods and return combined results.
        
        Args:
            peptides (list): List of peptides to analyze
            allele (str): MHC allele
            length (int): Peptide length
            
        Returns:
            pd.DataFrame: DataFrame with combined analysis results
        """
        logger.info(f"Analyzing {len(peptides)} peptides with allele {allele} using both EL and BA methods")
        
        # First get EL (eluted ligand) predictions
        el_results = self.predict_binding_api(peptides, allele, length, method="netmhcpan_el")
        
        # Then get BA (binding affinity) predictions
        ba_results = self.predict_binding_api(peptides, allele, length, method="netmhcpan_ba")
        
        if el_results.empty or ba_results.empty:
            logger.error("Failed to get predictions from one or both methods")
            if not el_results.empty:
                return el_results
            elif not ba_results.empty:
                return ba_results
            else:
                return pd.DataFrame()
        
        # Merge results
        # Keep only the necessary columns from each method
        try:
            # Format the allele name for immunogenicity calculation
            formatted_allele = allele.replace("*", "").replace(":", "")
            
            # Select columns from el_results
            if 'netmhcpan_el_score' in el_results.columns:
                el_subset = el_results[['peptide', 'netmhcpan_el_score', 'netmhcpan_el_rank']]
                el_subset = el_subset.rename(columns={
                    'netmhcpan_el_score': 'score',
                    'netmhcpan_el_rank': 'percentile_rank'
                })
            else:
                logger.warning("netmhcpan_el_score not found in EL results")
                el_subset = el_results[['peptide']]
                if 'score' in el_results.columns:
                    el_subset['score'] = el_results['score']
                if 'percentile_rank' in el_results.columns:
                    el_subset['percentile_rank'] = el_results['percentile_rank']
                
            # Select columns from ba_results
            if 'netmhcpan_ba_ic50' in ba_results.columns:
                ba_subset = ba_results[['peptide', 'netmhcpan_ba_ic50']]
                ba_subset = ba_subset.rename(columns={'netmhcpan_ba_ic50': 'ic50'})
            else:
                logger.warning("netmhcpan_ba_ic50 not found in BA results")
                ba_subset = ba_results[['peptide']]
                if 'ic50' in ba_results.columns:
                    ba_subset['ic50'] = ba_results['ic50']
                elif 'netmhcpan_ic50' in ba_results.columns:
                    ba_subset['ic50'] = ba_results['netmhcpan_ic50']
                
            # Merge the two dataframes on peptide
            combined = pd.merge(el_subset, ba_subset, on='peptide', how='outer')
            
            # Add allele information
            combined['allele'] = allele
            
            # Add immunogenicity scores
            combined = self.add_immunogenicity(combined, formatted_allele)
            
            # Fill any missing values in score/percentile_rank/ic50
            if 'score' not in combined.columns:
                combined['score'] = np.nan
            if 'percentile_rank' not in combined.columns:
                combined['percentile_rank'] = np.nan
            if 'ic50' not in combined.columns:
                combined['ic50'] = np.nan
                
            # Select and order final columns
            final_columns = ['peptide', 'allele', 'score', 'percentile_rank', 'immunogenicity', 'ic50']
            final_df = combined[final_columns]
            
            # Save combined results
            combined_csv = os.path.join(self.output_dir, f"combined_analysis_{allele.replace('*', '').replace(':', '')}.csv")
            self.save_to_csv(final_df, combined_csv)
            logger.info(f"Combined analysis saved to: {combined_csv}")
            
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
    
    def analyze_peptides_batch(self, peptides_list, alleles, length=9, batch_size=50):
        """
        Analyze batches of peptides across multiple alleles.
        
        Args:
            peptides_list (list): List of peptides
            alleles (list): List of MHC alleles
            length (int): Peptide length
            batch_size (int): Batch size
            
        Returns:
            dict: Dictionary of DataFrames with results for each allele
        """
        results = {}
        
        # Split peptides into batches to avoid overloading the API
        for i in range(0, len(peptides_list), batch_size):
            batch_peptides = peptides_list[i:i+batch_size]
            logger.info(f"Processing batch {i//batch_size + 1} with {len(batch_peptides)} peptides")
            
            # For each allele, analyze peptides
            for allele in alleles:
                logger.info(f"Analyzing peptides for allele: {allele}")
                batch_result = self.analyze_peptides(batch_peptides, allele, length)
                
                # Add results to dictionary
                if allele not in results:
                    results[allele] = batch_result
                else:
                    results[allele] = pd.concat([results[allele], batch_result], ignore_index=True)
        
        # Save combined results for each allele
        for allele, df in results.items():
            if not df.empty:
                combined_csv = os.path.join(self.output_dir, f"all_batches_{allele.replace('*', '').replace(':', '')}.csv")
                self.save_to_csv(df, combined_csv)
                logger.info(f"All batch results saved to: {combined_csv}")
        
        return results
    
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
            logger.info(f"Identified {len(binders)} binders that meet all thresholds")
            return binders
        else:
            logger.warning("Could not filter binders: no appropriate columns for filtering")
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
