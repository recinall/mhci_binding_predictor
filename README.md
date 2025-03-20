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

## Specifiche delle classi

### Immunogenicity

La classe `Immunogenicity` implementa l'algoritmo di predizione dell'immunogenicità descritto in Calis et al. (2013).

**Attributi principali:**
- `allele_dict`: Dizionario degli alleli disponibili con le posizioni mascherate per ciascuno
- `immunoscale`: Scala di immunogenicità per gli amminoacidi (valori derivati da Calis et al.)
- `immunoweight`: Pesi per le diverse posizioni nel peptide

**Metodi principali:**
- `predict_immunogenicity(peptides, custom_mask=None, allele=None)`: Predice l'immunogenicità di una lista di peptidi
- `validate_peptides(peptides, custom_mask=None, allele=None)`: Valida i peptidi prima della predizione
- `save_results_to_csv(results, filename)`: Salva i risultati in un file CSV
- `get_available_alleles()`: Restituisce la lista degli alleli disponibili
- `print_available_alleles()`: Stampa la lista degli alleli disponibili

**Funzionamento:**
La classe calcola l'immunogenicità di un peptide assegnando punteggi a ciascun amminoacido in base alla sua posizione nel peptide, escludendo le posizioni mascherate che sono tipicamente coinvolte nel binding MHC. I pesi delle posizioni e i valori degli amminoacidi sono derivati dall'analisi statistica di Calis et al. (2013).

### PeptideGenerator

La classe `PeptideGenerator` si occupa della generazione di peptidi 9-mer casuali da sequenze proteiche.

**Funzioni principali:**
- `download_swissprot_data()`: Scarica i dati da SwissProt o genera sequenze proteiche sintetiche
- `generate_9mers(protein_sequences, num_peptides)`: Genera peptidi 9-mer casuali dalle sequenze proteiche
- `save_to_csv(peptides, filename)`: Salva i peptidi generati in un file CSV
- `generate_peptides_pipeline(num_peptides, output_file)`: Pipeline completa per la generazione di peptidi

**Funzionamento:**
Il modulo scarica sequenze proteiche da SwissProt (o genera sequenze sintetiche come fallback) e poi estrae casualmente frammenti di 9 amminoacidi da queste sequenze, garantendo che contengano solo amminoacidi standard.

### PeptideAnalyzer

Il modulo `PeptideAnalyzer` gestisce l'analisi di binding MHC-I tramite il servizio IEDB.

**Funzioni principali:**
- `send_iedb_request(sequence_text, allele)`: Invia una richiesta al servizio IEDB
- `parse_iedb_response(response_text)`: Analizza la risposta del servizio IEDB
- `analyze_peptides(peptides, allele_list, batch_size, output_csv)`: Analizza una lista di peptidi
- `process_peptides_from_csv(input_csv, output_csv, allele_list, batch_size)`: Analizza peptidi da un file CSV

**Funzionamento:**
Il modulo invia richieste al servizio web IEDB per predire il binding dei peptidi agli alleli MHC-I specificati. Processa i risultati e li restituisce in un formato strutturato, includendo lo score di binding e il percentile rank per ciascun peptide.

### SequenceVariants

Il modulo `SequenceVariants` si occupa della generazione di varianti di sequenze peptidiche.

**Funzioni principali:**
- `generate_sequence_variants(sequenza)`: Genera tutte le possibili varianti di una sequenza
- `save_variants_to_csv(peptides, filename)`: Salva le varianti generate in un file CSV
- `generate_all_variants(sequences, output_dir, combined_file)`: Genera tutte le varianti per un insieme di sequenze

**Funzionamento:**
Il modulo prende in input una sequenza rappresentata come lista di liste, dove ogni lista interna contiene le possibili varianti per quella posizione, e genera tutte le possibili combinazioni di amminoacidi.

### Visualization

Il modulo `Visualization` fornisce funzioni per la visualizzazione dei risultati.

**Funzioni principali:**
- `plot_score_distribution(results, output_file)`: Grafico della distribuzione degli score
- `plot_category_distribution(results, output_file)`: Grafico della distribuzione delle categorie
- `plot_immunogenicity_correlation(results, output_file)`: Grafico della correlazione tra immunogenicità e binding
- `plot_percentile_distribution(results, output_file)`: Grafico della distribuzione dei percentile rank
- `plot_allele_comparison(results, output_file)`: Grafico di confronto tra alleli

**Funzionamento:**
Il modulo utilizza matplotlib e seaborn per creare visualizzazioni grafiche dei risultati dell'analisi, facilitando l'interpretazione dei dati e l'identificazione di pattern.

### CombinedResult

La classe `CombinedResult` gestisce i risultati combinati di binding MHC-I e immunogenicità.

**Attributi principali:**
- `peptide`: Sequenza peptidica
- `allele`: Allele HLA
- `binding_score`: Score di binding MHC-I
- `percentile_rank`: Percentile rank del binding MHC-I
- `immunogenicity_score`: Score di immunogenicità
- `composite_score`: Punteggio composito calcolato
- `category`: Categoria assegnata al peptide

**Metodi principali:**
- `_calculate_composite_score()`: Calcola il punteggio composito
- `_determine_category()`: Determina la categoria del peptide
- `to_dict()`: Converte l'oggetto in un dizionario
- `from_dict(data)`: Crea un oggetto CombinedResult da un dizionario
- `from_results_list(results_list)`: Crea una lista di oggetti CombinedResult da una lista di dizionari
- `add_immunogenicity_scores(results, output_csv, filtered_csv, ranked_csv)`: Aggiunge il punteggio di immunogenicità ai risultati
- `filter_by_percentile(results, threshold, operator)`: Filtra i risultati in base al percentile rank
- `filter_by_immunogenicity(results, threshold, operator)`: Filtra i risultati in base all'immunogenicità
- `filter_combined(results, percentile_threshold, percentile_operator, immunogenicity_threshold, immunogenicity_operator)`: Filtra i risultati con criteri combinati

**Funzionamento:**
La classe integra i risultati dell'analisi di binding MHC-I con i punteggi di immunogenicità, calcolando un punteggio composito e assegnando una categoria a ciascun peptide in base a criteri predefiniti.

## Calcolo del punteggio composito

Il punteggio composito combina il binding MHC-I e l'immunogenicità per fornire una valutazione complessiva del potenziale di un peptide come epitopo.

### Formula del punteggio composito

```
Punteggio composito = (immunogenicity_score * 0.5) + ((1 - percentile_rank/100) * 0.3) + (binding_score * 0.2)
```

Questa formula è stata progettata per:
1. Dare maggior peso all'immunogenicità (50%)
2. Valorizzare i peptidi con alto binding (percentile_rank basso) (30%)
3. Considerare anche il punteggio di binding assoluto (20%)
4. Produrre un valore facilmente interpretabile

### Componenti della formula:
- **Componente di immunogenicità**: `(immunogenicity_score * 0.5)`
  - Contribuisce per il 50% al punteggio finale
  - Varia da -0.5 a 0.5 (poiché immunogenicity_score varia da -1 a +1)

- **Componente di percentile rank**: `((1 - percentile_rank/100) * 0.3)`
  - Contribuisce per il 30% al punteggio finale
  - Varia da 0 a 0.3
  - Vale 0.3 quando percentile_rank = 0 (binding ottimale)
  - Vale 0 quando percentile_rank = 100 (binding pessimo)

- **Componente di binding score**: `(binding_score * 0.2)`
  - Contribuisce per il 20% al punteggio finale
  - Varia da 0 a 0.2 (poiché binding_score varia da 0 a 1)

### Intervallo di valori:
- Il punteggio composito varia da -0.5 a 1:
  - -0.5: peptide con immunogenicità estremamente negativa (score = -1), indipendentemente dal binding
  - 0: peptide con immunogenicità neutra (score = 0), pessimo binding (percentile_rank = 100) e binding score = 0
  - 1: peptide con massima immunogenicità positiva (score = 1), binding perfetto (percentile_rank = 0) e binding score = 1

### Sistema di categorizzazione

In base al punteggio composito e ad altri parametri specifici, i peptidi vengono categorizzati come:

| Categoria | Criteri |
|-----------|---------|
| **Eccellente** | percentile_rank < 0.1 E immunogenicity_score > 0.3 E binding_score > 0.95 |
| **Buono** | percentile_rank < 0.5 E immunogenicity_score > 0 E binding_score > 0.9 |
| **Da considerare** | percentile_rank < 1.0 E immunogenicity_score > 0 E binding_score > 0.8 |
| **Da scartare** | Tutti gli altri peptidi |

### Implementazione

La classe `CombinedResult` gestisce automaticamente il calcolo del punteggio composito e la categorizzazione dei peptidi, fornendo un'interfaccia unificata per l'analisi combinata di binding MHC-I e immunogenicità.

Il metodo `_calculate_composite_score()` implementa la formula sopra descritta, mentre il metodo `_determine_category()` assegna la categoria appropriata in base ai criteri definiti.

### Esempio di calcolo:

Per un peptide con:
- percentile_rank = 0.05 (eccellente binding)
- immunogenicity_score = 0.4 (buona immunogenicità)
- binding_score = 0.98 (alto binding)

Il punteggio composito sarà:
```
(0.4 * 0.5) + ((1 - 0.05/100) * 0.3) + (0.98 * 0.2) = 0.2 + 0.2999 + 0.196 = 0.6959
```

Questo punteggio, insieme agli altri parametri, classificherebbe il peptide nella categoria "Eccellente".

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
from peptide_analysis import load_results_from_csv, filter_results_by_percentile, filter_results_by_immunogenicity, filter_results_combined

# Carica i risultati da un file CSV
results = load_results_from_csv("results.csv")

# Filtra i risultati con percentile rank < 1.0
filtered_by_percentile = filter_results_by_percentile(results, threshold=1.0, operator="<")

# Filtra i risultati con immunogenicity score > 0
filtered_by_immunogenicity = filter_results_by_immunogenicity(results, threshold=0.0, operator=">")

# Filtra i risultati con entrambi i criteri
filtered_combined = filter_results_combined(
    results,
    percentile_threshold=0.5,
    percentile_operator="<",
    immunogenicity_threshold=0.0,
    immunogenicity_operator=">"
)
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

# Aggiungere l'immunogenicità a risultati esistenti
from peptide_analysis import load_results_from_csv, add_immunogenicity_to_results

# Carica i risultati da un file CSV
results = load_results_from_csv("results.csv")

# Aggiungi l'immunogenicità ai risultati
results_with_immuno, filtered_results, ranked_results = add_immunogenicity_to_results(
    results, 
    output_dir="immunogenicity_output"
)

# Utilizzo della classe CombinedResult per gestire i risultati combinati
from peptide_analysis import CombinedResult

# Crea un oggetto CombinedResult
combined_result = CombinedResult(
    peptide="GILGFVFTL",
    allele="HLA-A*02:01",
    binding_score=0.98,
    percentile_rank=0.05,
    immunogenicity_score=0.45
)

# Ottieni il punteggio composito e la categoria
print(f"Punteggio composito: {combined_result.composite_score}")
print(f"Categoria: {combined_result.category}")

# Converti in dizionario
result_dict = combined_result.to_dict()
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

# Filtraggio combinato (percentile rank e immunogenicità)
python peptide_analysis/examples/filter_combined_results.py --results-csv output/complete_results.csv --percentile-threshold 0.5 --percentile-operator "<" --immunogenicity-threshold 0.0 --immunogenicity-operator ">" --output-dir filtered_combined

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
- **immunogenicity.py**: Predizione dell'immunogenicità dei peptidi
- **combined_result.py**: Gestione dei risultati combinati di binding MHC-I e immunogenicità

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
