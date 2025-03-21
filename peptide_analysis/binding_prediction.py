"""
Modulo per la predizione del binding MHC-I utilizzando l'implementazione di IEDB.

Questo modulo fornisce funzioni per predire il binding di peptidi a molecole MHC-I
utilizzando l'algoritmo di predizione implementato in predict_binding.py.
"""

import os
import sys
import tempfile
import subprocess
import pandas as pd
import logging

# Configurazione del logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class BindingPredictor:
    """
    Classe per la predizione del binding MHC-I utilizzando l'implementazione di IEDB.
    
    Questa classe fornisce metodi per predire il binding di peptidi a molecole MHC-I
    utilizzando l'algoritmo di predizione implementato in predict_binding.py.
    """
    
    def __init__(self, predict_binding_path=None):
        """
        Inizializza il predittore di binding.
        
        Args:
            predict_binding_path: Percorso al file predict_binding.py (opzionale)
        """
        self.predict_binding_path = predict_binding_path
        if not self.predict_binding_path:
            # Cerca il file predict_binding.py nella directory corrente o in quella superiore
            current_dir = os.path.dirname(os.path.abspath(__file__))
            possible_paths = [
                os.path.join(current_dir, "predict_binding.py"),
                os.path.join(current_dir, "..", "predict_binding.py")
            ]
            
            for path in possible_paths:
                if os.path.exists(path):
                    self.predict_binding_path = path
                    break
        
        if not self.predict_binding_path or not os.path.exists(self.predict_binding_path):
            logger.warning("File predict_binding.py non trovato. Specificare il percorso corretto.")
    
    def predict_binding(self, peptides, allele="HLA-A*02:01", method="netmhcpan_el"):
        """
        Predice il binding di peptidi a una molecola MHC-I.
        
        Args:
            peptides: Lista di peptidi o singolo peptide
            allele: Allele MHC-I (default: HLA-A*02:01)
            method: Metodo di predizione (default: netmhcpan_el)
            
        Returns:
            DataFrame con i risultati della predizione
        """
        if not self.predict_binding_path:
            raise FileNotFoundError("File predict_binding.py non trovato. Specificare il percorso corretto.")
        
        if isinstance(peptides, str):
            peptides = [peptides]
        
        # Verifica che tutti i peptidi abbiano la stessa lunghezza
        peptide_length = len(peptides[0])
        if not all(len(p) == peptide_length for p in peptides):
            raise ValueError("Tutti i peptidi devono avere la stessa lunghezza")
        
        # Crea un file temporaneo con i peptidi in formato FASTA
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            temp_file.write(">Peptides\n")
            for peptide in peptides:
                temp_file.write(f"{peptide}\n")
            temp_file_path = temp_file.name
        
        try:
            # Esegui predict_binding.py
            cmd = [
                sys.executable,
                self.predict_binding_path,
                method,
                allele,
                str(peptide_length),
                temp_file_path
            ]
            
            logger.debug(f"Esecuzione comando: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Analizza l'output
            output_lines = result.stdout.strip().split('\n')
            
            # La prima riga contiene le intestazioni
            headers = output_lines[0].split('\t')
            
            # Le righe successive contengono i dati
            data = []
            for line in output_lines[1:]:
                if line.startswith('#'):  # Ignora le righe di commento
                    continue
                data.append(line.split('\t'))
            
            # Crea un DataFrame
            df = pd.DataFrame(data, columns=headers)
            
            # Converti le colonne numeriche
            numeric_columns = ['start', 'end', 'length']
            for col in numeric_columns:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col])
            
            # Converti le colonne di punteggio
            score_columns = [col for col in df.columns if 'score' in col.lower() or 'rank' in col.lower() or 'ic50' in col.lower()]
            for col in score_columns:
                if col in df.columns and df[col].dtype == object:
                    # Gestisci i valori '-' convertendoli in NaN
                    df[col] = df[col].replace('-', float('nan'))
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            
            return df
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Errore nell'esecuzione di predict_binding.py: {e}")
            logger.error(f"Output di errore: {e.stderr}")
            raise RuntimeError(f"Errore nella predizione del binding: {e}")
        
        finally:
            # Rimuovi il file temporaneo
            if os.path.exists(temp_file_path):
                os.unlink(temp_file_path)

def predict_binding(peptides, allele="HLA-A*02:01", method="netmhcpan_el", predict_binding_path=None):
    """
    Funzione di utilità per predire il binding di peptidi a una molecola MHC-I.
    
    Args:
        peptides: Lista di peptidi o singolo peptide
        allele: Allele MHC-I (default: HLA-A*02:01)
        method: Metodo di predizione (default: netmhcpan_el)
        predict_binding_path: Percorso al file predict_binding.py (opzionale)
        
    Returns:
        DataFrame con i risultati della predizione
    """
    predictor = BindingPredictor(predict_binding_path)
    return predictor.predict_binding(peptides, allele, method)

def get_available_methods(predict_binding_path=None):
    """
    Ottiene i metodi di predizione disponibili.
    
    Args:
        predict_binding_path: Percorso al file predict_binding.py (opzionale)
        
    Returns:
        Lista dei metodi disponibili
    """
    predictor = BindingPredictor(predict_binding_path)
    if not predictor.predict_binding_path:
        raise FileNotFoundError("File predict_binding.py non trovato. Specificare il percorso corretto.")
    
    cmd = [
        sys.executable,
        predictor.predict_binding_path,
        "method"
    ]
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )
    
    output_lines = result.stdout.strip().split('\n')
    methods = []
    capture = False
    
    for line in output_lines:
        if line == "MHC-I prediction methods:":
            capture = True
            continue
        elif line == "-------------------------":
            continue
        elif not line.strip():
            capture = False
        elif capture:
            methods.append(line.strip())
    
    return methods

def get_available_alleles(method="netmhcpan_el", predict_binding_path=None):
    """
    Ottiene gli alleli disponibili per un metodo specifico.
    
    Args:
        method: Metodo di predizione (default: netmhcpan_el)
        predict_binding_path: Percorso al file predict_binding.py (opzionale)
        
    Returns:
        DataFrame con gli alleli disponibili
    """
    predictor = BindingPredictor(predict_binding_path)
    if not predictor.predict_binding_path:
        raise FileNotFoundError("File predict_binding.py non trovato. Specificare il percorso corretto.")
    
    cmd = [
        sys.executable,
        predictor.predict_binding_path,
        method,
        "mhc"
    ]
    
    try:
        # Per i metodi con molti alleli, dobbiamo rispondere "yes" alla domanda
        import pexpect
        
        child = pexpect.spawn(' '.join(cmd))
        i = child.expect(['Do you still want to print it', pexpect.EOF, pexpect.TIMEOUT], timeout=5)
        
        if i == 0:
            child.sendline('y')
            child.expect(pexpect.EOF)
            output = child.before.decode('utf-8')
        else:
            # Se non c'è la domanda, eseguiamo normalmente
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            output = result.stdout
    except ImportError:
        # Se pexpect non è disponibile, proviamo con subprocess
        result = subprocess.run(
            cmd,
            input="y\n",
            capture_output=True,
            text=True
        )
        output = result.stdout
    
    output_lines = output.strip().split('\n')
    
    # Trova l'indice della riga di intestazione
    header_index = -1
    for i, line in enumerate(output_lines):
        if "Species" in line and "MHC" in line and "PeptideLength" in line:
            header_index = i
            break
    
    if header_index == -1:
        return pd.DataFrame()
    
    # Estrai i dati
    data = []
    for line in output_lines[header_index+1:]:
        if line.strip():
            parts = line.split('\t')
            if len(parts) >= 3:
                species = parts[0].strip()
                mhc = parts[1].strip()
                length = parts[2].strip()
                data.append([species, mhc, length])
    
    return pd.DataFrame(data, columns=["Species", "MHC", "PeptideLength"])
