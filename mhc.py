#!/usr/bin/env python3
"""
Optimized MHC-I Binding Predictor
Unified tool for peptide-MHC binding prediction using IEDB API
"""

import os
import logging
import tempfile
import numpy as np
import pandas as pd
import requests
import time
import click
from typing import List, Union, Optional, Dict, Any
import re
import textwrap

# Logger configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("IEDBBindingPredictor")

class IEDBBindingPredictor:
    """
    Optimized class for predicting peptide binding with MHC alleles using the IEDB API.
    """
    
    def __init__(self, output_dir="./output", csv_separator=",", decimal_separator="."):
        """Initialize the binding predictor."""
        self.output_dir = output_dir
        self.csv_separator = csv_separator
        self.decimal_separator = decimal_separator
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Immunogenicity calculation variables
        self.immunoscale = {
            "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110, 
            "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036, "M": -0.570, "N": -0.021, 
            "P": -0.036, "Q": -0.376, "R": 0.168, "S": -0.537, "T": 0.126, "V": 0.134, 
            "W": 0.719, "Y": -0.012
        }
        self.immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]
        
        # Allele anchor positions
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
        
        logger.info(f"Initialized binding predictor")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"CSV separator: '{self.csv_separator}', Decimal separator: '{self.decimal_separator}'")
    
    def _normalize_column_names(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Normalize column names from IEDB API response to standardized names.
        """
        if df.empty:
            return df
        
        # Log current columns for debugging
        logger.info(f"API response columns: {list(df.columns)}")
        
        # Common column name mappings from IEDB API
        column_mappings = {
            # Peptide sequence columns
            'seq': 'peptide',
            'sequence': 'peptide',
            'epitope': 'peptide',
            'peptide_seq': 'peptide',
            
            # Allele columns
            'mhc': 'allele',
            'hla': 'allele',
            'mhc_allele': 'allele',
            'hla_allele': 'allele',
            
            # Score columns
            'netmhcpan_el_score': 'el_score',
            'netmhcpan_ba_ic50': 'ic50',
            'netmhcpan_ba_score': 'ba_score',
            'score': 'el_score',  # Default score mapping
            
            # Percentile columns
            'percentile_rank': 'percentile_rank',
            'rank': 'percentile_rank',
            'el_rank': 'percentile_rank',
            
            # IC50 columns  
            'ic50': 'ic50',
            'ba_ic50': 'ic50'
        }
        
        # Create a copy to avoid modifying original
        df_normalized = df.copy()
        
        # Apply mappings (case-insensitive)
        for old_name, new_name in column_mappings.items():
            # Find matching column (case-insensitive)
            matching_cols = [col for col in df_normalized.columns if col.lower() == old_name.lower()]
            if matching_cols:
                df_normalized = df_normalized.rename(columns={matching_cols[0]: new_name})
        
        # If we still don't have peptide/allele columns, try to infer them
        if 'peptide' not in df_normalized.columns:
            # Look for any column that might contain peptide sequences
            potential_peptide_cols = [col for col in df_normalized.columns 
                                    if any(keyword in col.lower() 
                                          for keyword in ['seq', 'peptide', 'epitope'])]
            if potential_peptide_cols:
                df_normalized = df_normalized.rename(columns={potential_peptide_cols[0]: 'peptide'})
                logger.info(f"Mapped '{potential_peptide_cols[0]}' to 'peptide'")
        
        if 'allele' not in df_normalized.columns:
            # Look for any column that might contain allele names
            potential_allele_cols = [col for col in df_normalized.columns 
                                   if any(keyword in col.lower() 
                                         for keyword in ['allele', 'mhc', 'hla'])]
            if potential_allele_cols:
                df_normalized = df_normalized.rename(columns={potential_allele_cols[0]: 'allele'})
                logger.info(f"Mapped '{potential_allele_cols[0]}' to 'allele'")
        
        logger.info(f"Normalized columns: {list(df_normalized.columns)}")
        return df_normalized
    
    def _make_api_request(self, method: str, peptides: List[str], alleles: List[str], lengths: List[int]) -> pd.DataFrame:
        """
        Make a single optimized API request to IEDB.
        """
        url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
        
        # Format sequences as FASTA
        fasta_sequences = "\n".join([f">peptide{i+1}\n{p}" for i, p in enumerate(peptides)])
        
        # Prepare data for the request
        data = {
            "method": method,
            "sequence_text": fasta_sequences,
            "allele": ",".join(alleles),
            "length": ",".join(map(str, lengths))
        }
        
        logger.info(f"Making API request with method={method}, {len(peptides)} peptides, {len(alleles)} alleles")
        
        try:
            response = requests.post(url, data=data, timeout=60)
            
            if response.status_code != 200:
                logger.error(f"API request failed: {response.status_code} - {response.text}")
                return pd.DataFrame()
            
            # Parse response to DataFrame
            lines = response.text.strip().split("\n")
            if len(lines) < 2:
                logger.warning("No results from API")
                return pd.DataFrame()
            
            # Create temporary file to parse CSV properly
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.csv') as temp_file:
                temp_file.write(response.text)
                temp_file_path = temp_file.name
            
            try:
                df = pd.read_csv(temp_file_path, sep="\t")
                os.remove(temp_file_path)
                
                # Normalize column names
                df = self._normalize_column_names(df)
                
                return df
            except Exception as e:
                os.remove(temp_file_path)
                logger.error(f"Error parsing API response: {str(e)}")
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"API request error: {str(e)}")
            return pd.DataFrame()
    
    def calculate_immunogenicity_score(self, peptide: str, allele: str = None) -> float:
        """Calculate immunogenicity score for a peptide."""
        peptide = peptide.upper()
        peplen = len(peptide)
        cterm = peplen - 1
        score = 0
        
        # Determine mask positions
        if allele:
            clean_allele = allele.replace("*", "").replace(":", "")
            if clean_allele in self.allele_dict:
                mask_str = self.allele_dict[clean_allele].split(",")
                mask_positions = [int(x) - 1 for x in mask_str]  # Convert to 0-based
            else:
                mask_positions = [0, 1, cterm]  # Default mask
        else:
            mask_positions = [0, 1, cterm]  # Default mask
        
        # Adjust weights for longer peptides
        if peplen > 9:
            pepweight = self.immunoweight[:5] + ([0.30] * (peplen - 9)) + self.immunoweight[5:]
        else:
            pepweight = self.immunoweight[:peplen]
        
        try:
            for i, aa in enumerate(peptide):
                if aa not in self.immunoscale:
                    logger.warning(f"Invalid amino acid '{aa}' in peptide {peptide}")
                    return 0.0
                
                if i not in mask_positions and i < len(pepweight):
                    score += pepweight[i] * self.immunoscale[aa]
            
            return round(score, 5)
        
        except Exception as e:
            logger.error(f"Error calculating immunogenicity for {peptide}: {str(e)}")
            return 0.0
    
    def predict_comprehensive(self, peptides: List[str], alleles: List[str], lengths: List[int] = None, delay: float = 2.0) -> pd.DataFrame:
        """
        Make comprehensive predictions using both EL and BA methods with one request per allele.
        """
        if lengths is None:
            lengths = [9]
        
        # Convert single values to lists
        if isinstance(alleles, str):
            alleles = [alleles]
        if isinstance(lengths, int):
            lengths = [lengths]
        
        all_results = []
        
        logger.info(f"Processing {len(alleles)} alleles with {delay}s delay between requests...")
        
        # Process each allele separately to avoid API limits
        for i, allele in enumerate(alleles):
            logger.info(f"Processing allele {i+1}/{len(alleles)}: {allele}")
            
            # Add delay between requests (except for the first one)
            if i > 0:
                logger.info(f"Waiting {delay} seconds before next request...")
                time.sleep(delay)
            
            allele_results = []
            
            # Get EL predictions for this allele
            logger.info(f"Getting EL predictions for {allele}...")
            el_results = self._make_api_request("netmhcpan_el", peptides, [allele], lengths)
            
            if not el_results.empty:
                # Check if required columns exist
                if 'peptide' not in el_results.columns or 'allele' not in el_results.columns:
                    logger.warning(f"Missing required columns in EL results for {allele}. Available columns: {list(el_results.columns)}")
                    continue
                
                # Select relevant columns from EL results
                available_cols = ['allele', 'peptide', 'percentile_rank']
                if 'el_score' in el_results.columns:
                    available_cols.append('el_score')
                elif 'score' in el_results.columns:
                    available_cols.append('score')
                    el_results = el_results.rename(columns={'score': 'el_score'})
                    available_cols[-1] = 'el_score'
                
                el_data = el_results[[col for col in available_cols if col in el_results.columns]].copy()
                el_data['method'] = 'netmhcpan_el'
                
                # Add delay before BA request
                logger.info(f"Waiting {delay} seconds before BA request...")
                time.sleep(delay)
                
                # Get BA predictions for this allele
                logger.info(f"Getting BA predictions for {allele}...")
                ba_results = self._make_api_request("netmhcpan_ba", peptides, [allele], lengths)
                
                if not ba_results.empty:
                    # Check if required columns exist
                    if 'peptide' in ba_results.columns and 'allele' in ba_results.columns:
                        # Select relevant columns from BA results
                        ba_cols = ['allele', 'peptide']
                        if 'ic50' in ba_results.columns:
                            ba_cols.append('ic50')
                        
                        ba_data = ba_results[[col for col in ba_cols if col in ba_results.columns]].copy()
                        
                        # Merge BA results with EL results
                        if 'ic50' in ba_data.columns:
                            el_data = el_data.merge(
                                ba_data[['allele', 'peptide', 'ic50']], 
                                on=['allele', 'peptide'], 
                                how='left'
                            )
                
                all_results.append(el_data)
                logger.info(f"Completed predictions for {allele}: {len(el_data)} results")
            else:
                logger.warning(f"No EL results obtained for allele {allele}")
        
        # Combine all results
        if all_results:
            combined_df = pd.concat(all_results, ignore_index=True)
            logger.info(f"Combined results from all alleles: {len(combined_df)} total results")
        else:
            logger.error("No results obtained from API for any allele")
            return pd.DataFrame()
        
        # Verify we have the required columns before calculating immunogenicity
        if 'peptide' not in combined_df.columns or 'allele' not in combined_df.columns:
            logger.error(f"Missing required columns for immunogenicity calculation. Available: {list(combined_df.columns)}")
            return combined_df
        
        # Add immunogenicity scores
        logger.info("Calculating immunogenicity scores...")
        try:
            combined_df['immunogenicity'] = combined_df.apply(
                lambda row: self.calculate_immunogenicity_score(row['peptide'], row['allele']), 
                axis=1
            )
        except Exception as e:
            logger.error(f"Error calculating immunogenicity scores: {str(e)}")
            combined_df['immunogenicity'] = 0.0
        
        # Standardize column names and structure
        final_df = self._standardize_columns(combined_df)
        
        return final_df
    
    def _standardize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize column names and ensure all required columns exist."""
        if df.empty:
            return df
        
        # Ensure required columns exist with default values
        required_columns = ['peptide', 'allele', 'percentile_rank', 'el_score', 'ic50', 'immunogenicity']
        
        for col in required_columns:
            if col not in df.columns:
                df[col] = np.nan
        
        # Select final columns in order (only those that exist)
        final_columns = [col for col in ['peptide', 'allele', 'el_score', 'percentile_rank', 'ic50', 'immunogenicity'] 
                        if col in df.columns]
        
        return df[final_columns]
    
    def save_to_csv(self, df: pd.DataFrame, file_path: str):
        """Save DataFrame to CSV with configured separators."""
        if df.empty:
            logger.warning(f"Empty DataFrame, not saving to {file_path}")
            return
            
        try:
            df.to_csv(file_path, 
                     sep=self.csv_separator, 
                     decimal=self.decimal_separator, 
                     index=False)
            logger.info(f"Saved {len(df)} results to {file_path}")
        except Exception as e:
            logger.error(f"Error saving to CSV: {str(e)}")
    
    def filter_binders(self, df: pd.DataFrame, 
                      el_score_threshold: float = None,
                      percentile_threshold: float = None,
                      ic50_threshold: float = None,
                      immunogenicity_threshold: float = None) -> pd.DataFrame:
        """Filter binders based on multiple criteria."""
        if df.empty:
            return df
        
        # Convert numeric columns
        numeric_cols = ['el_score', 'percentile_rank', 'ic50', 'immunogenicity']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Apply filters
        mask = pd.Series([True] * len(df))
        
        if el_score_threshold is not None and 'el_score' in df.columns:
            mask &= (df['el_score'] >= el_score_threshold)
            
        if percentile_threshold is not None and 'percentile_rank' in df.columns:
            mask &= (df['percentile_rank'] <= percentile_threshold)
            
        if ic50_threshold is not None and 'ic50' in df.columns:
            mask &= (df['ic50'] <= ic50_threshold)
            
        if immunogenicity_threshold is not None and 'immunogenicity' in df.columns:
            mask &= (df['immunogenicity'] >= immunogenicity_threshold)
        
        filtered_df = df[mask].copy()
        logger.info(f"Filtered {len(df)} results to {len(filtered_df)} binders")
        
        return filtered_df
    
    def generate_variants(self, pattern: str) -> List[str]:
        """Generate peptide variants from a pattern like A[CD]E[FY]GH."""
        variants = []
        
        def expand_pattern(current: str, remaining: str):
            if not remaining:
                variants.append(current)
                return
            
            if remaining[0] == '[':
                # Find closing bracket
                close_idx = remaining.find(']')
                if close_idx == -1:
                    # Invalid pattern, treat as literal
                    expand_pattern(current + remaining[0], remaining[1:])
                    return
                
                # Extract options
                options = remaining[1:close_idx]
                for option in options:
                    expand_pattern(current + option, remaining[close_idx+1:])
            else:
                # Regular character
                expand_pattern(current + remaining[0], remaining[1:])
        
        expand_pattern("", pattern)
        
        if not variants:
            variants = [pattern]  # Fallback to original pattern
        
        logger.info(f"Generated {len(variants)} variants from pattern: {pattern}")
        return variants

# CLI Interface
@click.group(epilog=textwrap.dedent('''\
    ESEMPI DETTAGLIATI:
    
      PREDIZIONE:
        Predizione base con delay personalizzato:
        python mhci_predictor.py predict --peptides "SIINFEKL,RAKFKQLL" --alleles "HLA-A*02:01,HLA-B*07:02" --delay 3.0
    
        Predizione da file di peptidi:
        python mhci_predictor.py predict --peptides peptide_list.txt --alleles "HLA-A*02:01" --output risultati.csv
    
      ANALISI PATTERN:
        Analisi pattern con lunghezze multiple:
        python mhci_predictor.py analyze --pattern "A[CD]E[FY]GH" --alleles "HLA-A*02:01" --lengths "8,9,10" --output analisi.csv
    
      FILTRO:
        Filtro con soglie multiple:
        python mhci_predictor.py filter --input risultati.csv --percentile 2.0 --ic50 500 --immunogenicity 0.5
    
      VARIANTI:
        Generazione varianti con salvataggio su file:
        python mhci_predictor.py variants "A[CD]E[FY]GH" --output varianti.txt

    PARAMETRI GLOBALI:
      --output-dir:  Directory di output (default: ./output)
      --csv-sep:     Separatore colonne CSV (default: ',')
      --decimal-sep: Separatore decimali (default: '.')
'''))
@click.option('--output-dir', default='./output', help='Output directory')
@click.option('--csv-sep', default=',', help='CSV separator')
@click.option('--decimal-sep', default='.', help='Decimal separator')
@click.pass_context
def main(ctx, output_dir, csv_sep, decimal_sep):
    """Optimized MHC-I Binding Prediction Tool"""
    ctx.ensure_object(dict)
    ctx.obj['predictor'] = IEDBBindingPredictor(
        output_dir=output_dir,
        csv_separator=csv_sep,
        decimal_separator=decimal_sep
    )

@main.command()
@click.option('--peptides', required=True, help='''Percorso a file o lista di peptidi separati da virgola
              (formato file: un peptide per riga)''')
@click.option('--alleles', required=True, help='Lista di alleli MHC separati da virgola')
@click.option('--lengths', default='9', help='Lunghezze peptidi separate da virgola (default: 9)')
@click.option('--output', help='Percorso file output CSV')
@click.option('--delay', default=2.0, type=float, help='Delay tra richieste API in secondi (default: 2.0)')
@click.pass_context
def predict(ctx, peptides, alleles, lengths, output, delay):
    """Esegue predizioni di binding complete
    (metodi EL + BA + calcolo immunogenicita')
    """
    predictor = ctx.obj['predictor']
    
    # Parse peptides
    if os.path.exists(peptides):
        with open(peptides, 'r') as f:
            peptides_list = [line.strip() for line in f if line.strip()]
    else:
        peptides_list = [p.strip() for p in peptides.split(',')]
    
    # Parse alleles and lengths
    alleles_list = [a.strip() for a in alleles.split(',')]
    lengths_list = [int(l.strip()) for l in lengths.split(',')]
    
    # Run predictions
    results = predictor.predict_comprehensive(peptides_list, alleles_list, lengths_list, delay=delay)
    
    # Save results
    if output:
        predictor.save_to_csv(results, output)
    else:
        output_file = os.path.join(predictor.output_dir, "prediction_results.csv")
        predictor.save_to_csv(results, output_file)
    
    click.echo(f"✅ Predictions completed. Results saved to {output or output_file}")

@main.command()
@click.option('--pattern', required=True, help='Sequence pattern (e.g., A[CD]E[FY]GH)')
@click.option('--alleles', required=True, help='MHC alleles (comma-separated)')
@click.option('--lengths', default='9', help='Peptide lengths (comma-separated)')
@click.option('--output', required=True, help='Output CSV file')
@click.option('--delay', default=2.0, type=float, help='Delay between API requests (seconds)')
@click.pass_context
def analyze(ctx, pattern, alleles, lengths, output, delay):
    """Analyze peptide pattern variants"""
    predictor = ctx.obj['predictor']
    
    # Generate variants
    peptides_list = predictor.generate_variants(pattern)
    alleles_list = [a.strip() for a in alleles.split(',')]
    lengths_list = [int(l.strip()) for l in lengths.split(',')]
    
    # Run predictions
    results = predictor.predict_comprehensive(peptides_list, alleles_list, lengths_list, delay=delay)
    
    # Save results
    predictor.save_to_csv(results, output)
    click.echo(f"✅ Pattern analysis completed. {len(results)} results saved to {output}")

@main.command()
@click.option('--input', required=True, help='File CSV di input con risultati predizioni')
@click.option('--output', help='File output per risultati filtrati')
@click.option('--el-score', type=float, help='Soglia minima EL score')
@click.option('--percentile', type=float, help='Soglia massima percentile rank')
@click.option('--ic50', type=float, help='Soglia massima IC50 (nM)')
@click.option('--immunogenicity', type=float, help='Soglia minima immunogenicita')
@click.pass_context
def filter(ctx, input, output, el_score, percentile, ic50, immunogenicity):
    """Filtra risultati in base a soglie specificate
    (almeno una soglia deve essere fornita)
    """
    predictor = ctx.obj['predictor']
    
    # Load data
    df = pd.read_csv(input, sep=predictor.csv_separator, decimal=predictor.decimal_separator)
    
    # Apply filters
    filtered = predictor.filter_binders(
        df,
        el_score_threshold=el_score,
        percentile_threshold=percentile,
        ic50_threshold=ic50,
        immunogenicity_threshold=immunogenicity
    )
    
    # Output results
    if output:
        predictor.save_to_csv(filtered, output)
        click.echo(f"✅ Filtered {len(filtered)} results saved to {output}")
    else:
        click.echo(filtered.to_string(index=False))

@main.command()
@click.argument('pattern')
@click.option('--output', help='File output per le varianti')
def variants(pattern, output):
    """Genera tutte le varianti peptidiche da un pattern
    Pattern: sequenza con gruppi opzionali in parentesi quadre
    Esempio: "A[CD]E[FY]GH" genera ACEFGH, ACEYGH, ADFGH, ADYGH
    """
    predictor = IEDBBindingPredictor()
    variants_list = predictor.generate_variants(pattern)
    
    if output:
        with open(output, 'w') as f:
            f.write('\n'.join(variants_list))
        click.echo(f"✅ Generated {len(variants_list)} variants saved to {output}")
    else:
        for variant in variants_list:
            click.echo(variant)

if __name__ == '__main__':
    main()
