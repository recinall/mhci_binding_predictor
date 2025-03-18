"""
Modulo per la visualizzazione dei risultati dell'analisi dei peptidi.
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from io import StringIO

def plot_score_distribution(results, output_file=None):
    """
    Crea un grafico della distribuzione degli score per ogni allele.
    
    Parametri:
    results (list): Lista di risultati dell'analisi
    output_file (str): Nome del file di output (opzionale)
    
    Returns:
    matplotlib.figure.Figure: Figura creata
    """
    # Convertiamo i risultati in un DataFrame
    df = pd.DataFrame(results)
    
    # Creiamo la figura
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Creiamo un boxplot per ogni allele
    sns.boxplot(x='allele', y='score', data=df, ax=ax)
    
    # Aggiungiamo titolo e label
    ax.set_title('Distribuzione degli score per allele', fontsize=16)
    ax.set_xlabel('Allele', fontsize=14)
    ax.set_ylabel('Score', fontsize=14)
    
    # Ruotiamo le etichette dell'asse x
    plt.xticks(rotation=45, ha='right')
    
    # Aggiungiamo una griglia
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Aggiustiamo il layout
    plt.tight_layout()
    
    # Salviamo la figura se è stato specificato un file di output
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Grafico salvato in {output_file}")
    
    return fig

def plot_percentile_distribution(results, output_file=None):
    """
    Crea un grafico della distribuzione dei percentile rank per ogni allele.
    
    Parametri:
    results (list): Lista di risultati dell'analisi
    output_file (str): Nome del file di output (opzionale)
    
    Returns:
    matplotlib.figure.Figure: Figura creata
    """
    # Convertiamo i risultati in un DataFrame
    df = pd.DataFrame(results)
    
    # Creiamo la figura
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Creiamo un violin plot per ogni allele
    sns.violinplot(x='allele', y='percentile_rank', data=df, ax=ax)
    
    # Aggiungiamo titolo e label
    ax.set_title('Distribuzione dei percentile rank per allele', fontsize=16)
    ax.set_xlabel('Allele', fontsize=14)
    ax.set_ylabel('Percentile Rank', fontsize=14)
    
    # Ruotiamo le etichette dell'asse x
    plt.xticks(rotation=45, ha='right')
    
    # Aggiungiamo una griglia
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Aggiustiamo il layout
    plt.tight_layout()
    
    # Salviamo la figura se è stato specificato un file di output
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Grafico salvato in {output_file}")
    
    return fig

def plot_allele_comparison(results, output_file=None):
    """
    Crea un grafico di confronto tra gli alleli.
    
    Parametri:
    results (list): Lista di risultati dell'analisi
    output_file (str): Nome del file di output (opzionale)
    
    Returns:
    matplotlib.figure.Figure: Figura creata
    """
    # Convertiamo i risultati in un DataFrame
    df = pd.DataFrame(results)
    
    # Creiamo un pivot table per confrontare gli alleli
    pivot_df = df.pivot_table(
        index='peptide', 
        columns='allele', 
        values='percentile_rank',
        aggfunc='mean'
    ).reset_index()
    
    # Otteniamo gli alleli unici
    alleles = df['allele'].unique()
    
    # Creiamo la figura
    fig, axes = plt.subplots(1, len(alleles), figsize=(16, 6))
    
    # Se c'è un solo allele, convertiamo axes in una lista
    if len(alleles) == 1:
        axes = [axes]
    
    # Creiamo un istogramma per ogni allele
    for i, allele in enumerate(alleles):
        allele_data = df[df['allele'] == allele]['percentile_rank']
        sns.histplot(allele_data, ax=axes[i], kde=True)
        axes[i].set_title(f'Distribuzione per {allele}')
        axes[i].set_xlabel('Percentile Rank')
        axes[i].set_ylabel('Frequenza')
    
    # Aggiustiamo il layout
    plt.tight_layout()
    
    # Salviamo la figura se è stato specificato un file di output
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Grafico salvato in {output_file}")
    
    # Creiamo anche un grafico di correlazione se ci sono più alleli
    if len(alleles) > 1:
        # Calcoliamo la matrice di correlazione
        corr_matrix = pivot_df.drop('peptide', axis=1).corr()
        
        # Creiamo la figura
        fig2, ax2 = plt.subplots(figsize=(10, 8))
        
        # Creiamo una heatmap
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', ax=ax2)
        
        # Aggiungiamo titolo
        ax2.set_title('Correlazione tra alleli', fontsize=16)
        
        # Aggiustiamo il layout
        plt.tight_layout()
        
        # Salviamo la figura se è stato specificato un file di output
        if output_file:
            corr_output = output_file.replace('.', '_corr.')
            plt.savefig(corr_output, dpi=300, bbox_inches='tight')
            print(f"Grafico di correlazione salvato in {corr_output}")
    
    return fig
