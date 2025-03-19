# Peptide Analysis Library

Una libreria Python completa per la generazione, analisi e visualizzazione di peptidi 9-mer, con particolare attenzione all'analisi di binding MHC-I tramite il servizio IEDB e alla predizione dell'immunogenicità.

## Caratteristiche

- **Generazione di peptidi**: Creazione di peptidi 9-mer casuali da sequenze proteiche
- **Analisi di binding MHC-I**: Integrazione con il servizio IEDB per la predizione di binding
- **Predizione dell'immunogenicità**: Implementazione dell'algoritmo di Calis et al. (2013) per la predizione dell'immunogenicità dei peptidi
- **Categorizzazione degli epitopi**: Classificazione automatica dei peptidi in base a score, percentile rank e immunogenicità
- **Punteggio composito**: Calcolo di un punteggio che combina binding MHC-I e immunogenicità
- **Generazione di varianti di sequenze**: Creazione di tutte le possibili varianti di sequenze peptidiche
- **Visualizzazione dei risultati**: Grafici per l'analisi dei risultati (distribuzione degli score, percentile rank, ecc.)
- **Filtraggio dei risultati**: Possibilità di filtrare i risultati in base al percentile rank

## Installazione

```bash
# Clona il repository
git clone https://github.com/username/peptide_analysis.git
cd peptide_analysis

# Installa la libreria in modalità sviluppo
pip install -e .
```

## Utilizzo

### Generazione e analisi di peptidi casuali

```python
from peptide_analysis import generate_and_analyze_peptides

# Genera 1000 peptidi casuali e li analizza con gli alleli specificati
peptides, results = generate_and_analyze_peptides(
    num_peptides=1000,
    allele_list="HLA-A*01:01,HLA-A*02:01",
    batch_size=10,
    output_dir="output"
)
```

### Analisi completa con visualizzazione

```python
from peptide_analysis import run_complete_analysis

# Esegue un'analisi completa con visualizzazione
report = run_complete_analysis(
    input_csv="peptides.csv",  # Opzionale, se non specificato genera peptidi casuali
    num_peptides=1000,
    allele_list="HLA-A*01:01,HLA-A*02:01",
    batch_size=10,
    output_dir="output",
    percentile_threshold=10.0,
    percentile_operator="<="
)
```

### Generazione e analisi di varianti di sequenze

```python
from peptide_analysis import analyze_sequence_variants

# Definisci le sequenze con varianti
sequences = [
    [["S"], ["D","G","N"], ["P"], ["A","K","V"], ["R","C"], ["Y","A","H","N"], ["E","P","H"], ["F","H"], ["L"]],
    [["F","H"], ["L"], ["W","I","N"], ["G","V","L"], ["P","T","H"], ["R","N"], ["A","T"], ["L","V","H"], ["A","V","I"]]
]

# Genera e analizza tutte le varianti
report = analyze_sequence_variants(
    sequences=sequences,
    allele_list="HLA-A*01:01,HLA-A*02:01",
    batch_size=10,
    output_dir="variants_output",
    percentile_threshold=5.0,
    percentile_operator="<"
)
```

### Filtraggio dei risultati

```python
from peptide_analysis import load_results_from_csv, filter_results_by_percentile

# Carica i risultati da un file CSV
results = load_results_from_csv("results.csv")

# Filtra i risultati con percentile rank < 1.0
filtered_results = filter_results_by_percentile(results, threshold=1.0, operator="<")
```

### Predizione dell'immunogenicità e categorizzazione

```python
from peptide_analysis import predict_peptide_immunogenicity, print_available_immunogenicity_alleles

# Mostra gli alleli disponibili per la predizione dell'immunogenicità
print_available_immunogenicity_alleles()

# Predici l'immunogenicità di una lista di peptidi
peptides = ["GILGFVFTL", "NLVPMVATV", "CINGVCWTV"]
results = predict_peptide_immunogenicity(
    peptides=peptides,
    allele="HLA-A0201",  # Opzionale
    custom_mask="1,2,9",  # Opzionale
    output_csv="immunogenicity_results.csv"  # Opzionale
)

# Analisi completa dell'immunogenicità
from peptide_analysis import analyze_peptide_immunogenicity

report = analyze_peptide_immunogenicity(
    peptides=peptides,
    allele="HLA-A0201",
    output_dir="immunogenicity_output"
)

# Analisi completa con categorizzazione automatica
from peptide_analysis import run_complete_analysis

report = run_complete_analysis(
    peptides=peptides,
    allele_list="HLA-A*02:01",
    output_dir="complete_analysis",
    include_immunogenicity=True  # Abilita la categorizzazione automatica
)

# I risultati includeranno:
# - Punteggio di immunogenicità specifico per allele
# - Categoria (Eccellente, Buono, Da considerare, Da scartare)
# - Punteggio composito che combina binding e immunogenicità
```

### Visualizzazione dei risultati

```python
from peptide_analysis import (
    plot_score_distribution,
    plot_percentile_distribution,
    plot_allele_comparison
)

# Crea grafici per i risultati
plot_score_distribution(results, "score_distribution.png")
plot_percentile_distribution(results, "percentile_distribution.png")
plot_allele_comparison(results, "allele_comparison.png")
```

## Utilizzo da linea di comando

La libreria include anche script di esempio per l'utilizzo da linea di comando:

```bash
# Analisi completa di peptidi
python example.py --mode analyze --num-peptides 1000 --output-dir output

# Analisi di peptidi da un file CSV esistente
python example.py --mode analyze --input-csv peptides.csv --output-dir output

# Visualizzazione dei risultati di un'analisi precedente
python example.py --mode visualize --results-csv output/results.csv --output-dir output

# Generazione e analisi di varianti di sequenze
python peptide_analysis/examples/sequence_variants_example.py --output-dir variants_output

# Filtraggio dei risultati
python peptide_analysis/examples/filter_results.py --results-csv output/results.csv --percentile-threshold 1.0 --percentile-operator "<" --output-dir filtered_results

# Predizione dell'immunogenicità
python peptide_analysis/examples/immunogenicity_example.py --peptides GILGFVFTL NLVPMVATV CINGVCWTV --allele HLA-A0201 --output-dir immunogenicity_test

# Visualizzazione degli alleli disponibili per l'immunogenicità
python peptide_analysis/examples/immunogenicity_example.py --list-alleles
```

## Struttura della libreria

- **peptide_generator.py**: Generazione di peptidi casuali
- **peptide_analyzer.py**: Analisi di binding MHC-I tramite IEDB
- **sequence_variants.py**: Generazione di varianti di sequenze
- **visualization.py**: Visualizzazione dei risultati
- **utils.py**: Funzioni di utilità
- **main.py**: Funzioni principali che combinano le varie funzionalità

## Requisiti

- Python >= 3.6
- requests
- pandas
- matplotlib
- seaborn
- tqdm
- numpy
- openpyxl

## Licenza

MIT
