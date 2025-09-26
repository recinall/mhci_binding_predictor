# MHC-I Binding Predictor

Strumento completo per la predizione del legame peptide-MHC di classe I utilizzando l'API IEDB, con interfaccia grafica (GUI) e interfaccia a riga di comando (CLI).

## Panoramica

Questo progetto fornisce due modalità di utilizzo:
1. **GUI moderna** basata su PySide6 per un'esperienza utente intuitiva
2. **CLI potente** per automazione e elaborazione batch

Le predizioni vengono effettuate tramite l'API IEDB utilizzando i metodi NetMHCpan per:
- Score di legame EL (Elicitation Likelihood)
- Valori IC50 (Binding Affinity)
- Calcolo dell'immunogenicità
- Analisi di pattern peptidici

## Caratteristiche Principali

### Interfaccia Grafica (gui.py)
- **Predizioni dirette**: Inserimento di peptidi e alleli per ottenere score completi
- **Analisi pattern**: Generazione e analisi di varianti peptidiche da pattern (es. `A[CD]E[FY]GH`)
- **Filtri avanzati**: Filtraggio risultati basato su multiple soglie
- **Esportazione**: Salvataggio risultati in CSV con separatori personalizzabili
- **Interfaccia a schede**: Organizzazione chiara delle diverse funzionalità

### Interfaccia CLI (mhc.py)
- **Predizioni batch**: Elaborazione di liste di peptidi da file o input diretto
- **Analisi pattern**: Generazione varianti e predizioni combinate
- **Filtraggio risultati**: Applicazione soglie a file esistenti
- **Output flessibile**: Supporto per diversi separatori CSV e decimali

## Installazione

### Prerequisiti
- Python 3.7+
- Connessione internet (per l'API IEDB)

### Dipendenze
Installa le dipendenze richieste:

```bash
pip install pandas numpy requests PySide6 click
```

### Struttura del Progetto
```
.
├── gui.py          # Interfaccia grafica principale
├── mhc.py          # Motore di predizione e CLI
└── output/         # Directory per i risultati (creata automaticamente)
```

## Utilizzo della GUI

Avvia l'interfaccia grafica:

```bash
python gui.py
```

### Schede Disponibili

1. **Predictions** (Predizioni)
   - Inserisci peptidi (uno per riga) o carica da file
   - Specifica gli alleli MHC (separati da virgola)
   - Imposta lunghezze peptidiche e delay tra richieste
   - Clicca "Run Predictions" per avviare le predizioni

2. **Pattern Analysis** (Analisi Pattern)
   - Inserisci un pattern (es. `A[CD]E[FY]GH`)
   - Specifica alleli e parametri
   - Clicca "Analyze Pattern" per generare varianti ed eseguire predizioni

3. **Filter Results** (Filtra Risultati)
   - Applica filtri ai risultati correnti
   - Soglie disponibili: EL Score, Percentile Rank, IC50, Immunogenicity
   - Clicca "Apply Filters" per aggiornare la visualizzazione

### Funzionalità dei Risultati
- Tabella ordinabile con tutti i dati delle predizioni
- Salvataggio risultati completi o filtrati
- Esportazione con separatori personalizzabili
- Pulizia della visualizzazione

## Utilizzo della CLI

Il modulo `mhc.py` offre un'interfaccia a riga di comando completa:

### Comandi Disponibili

#### Predizioni
```bash
python mhc.py predict --peptides "SIINFEKL,RAKFKQLL" --alleles "HLA-A*02:01,HLA-B*07:02" --delay 3.0
```

```bash
python mhc.py predict --peptides peptide_list.txt --alleles "HLA-A*02:01" --output risultati.csv
```

#### Analisi Pattern
```bash
python mhc.py analyze --pattern "A[CD]E[FY]GH" --alleles "HLA-A*02:01" --lengths "8,9,10" --output analisi.csv
```

#### Filtraggio
```bash
python mhc.py filter --input risultati.csv --percentile 2.0 --ic50 500 --immunogenicity 0.5
```

#### Generazione Varianti
```bash
python mhc.py variants "A[CD]E[FY]GH" --output varianti.txt
```

### Opzioni Globali
- `--output-dir`: Directory di output (default: `./output`)
- `--csv-sep`: Separatore colonne CSV (default: `,`)
- `--decimal-sep`: Separatore decimali (default: `.`)

## Output delle Predizioni

I risultati includono:
- **peptide**: Sequenza peptidica
- **allele**: Allele MHC
- **el_score**: Score di legame EL (0-1, valori più alti indicano legame migliore)
- **percentile_rank**: Rank percentile (valori più bassi indicano legame migliore)
- **ic50**: Valore IC50 in nM (valori più bassi indicano legame migliore)
- **immunogenicity**: Score di immunogenicità calcolato

## Soglie Consigliate

Per identificare legami forti:
- **EL Score**: ≥ 0.5
- **Percentile Rank**: ≤ 2.0
- **IC50**: ≤ 500 nM
- **Immunogenicity**: ≥ 0.5 (score relativo)

## Note Tecniche

### API IEDB
- Le richieste vengono effettuate a `http://tools-cluster-interface.iedb.org/tools_api/mhci/`
- Vengono utilizzati i metodi `netmhcpan_el` e `netmhcpan_ba`
- È incluso un delay tra le richieste per evitare limitazioni API

### Immunogenicità
Lo score di immunogenicità viene calcolato basandosi su:
- Pesi aminoacidici specifici
- Posizioni di ancoraggio per allele MHC
- Mascheramento delle posizioni chiave

### Gestione Errori
- Il tool include logging completo per il debug
- Gestione degli errori di rete e API
- Validazione degli input

## Esempi Dettagliati

Esempi completi sono disponibili eseguendo:
```bash
python mhc.py --help
```

## Supporto

Per problemi o domande, verificare:
- La connessione internet
- Il formato degli input (peptidi e alleli)
- I log generati per errori specifici

## IEDB Data Attribution & License

This toolkit utilizes data from the Immune Epitope Database (IEDB). All IEDB data are released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license. By using this toolkit, you agree to:

1. **Credit IEDB**: Include attribution in publications or reports, e.g.:  
   *Data courtesy of the Immune Epitope Database (IEDB; https://www.iedb.org) under CC BY 4.0.*
2. **Link to License**: https://creativecommons.org/licenses/by/4.0/
3. **Indicate Modifications**: State if data have been altered from the original source.

---

## API Rate Limiting & Best Practices

To ensure responsible use of the public IEDB REST API, observe the following:

- **Throttle Requests**: Limit to a few requests per second per IP.  
- **Handle 429 Responses**: Implement exponential back-off on HTTP 429 (Too Many Requests) and respect any `Retry-After` headers.
- **Parallelism**: Use moderate concurrency (e.g. thread pool or async limited to ~5 parallel calls) to avoid server overload.

---

## Commercial Licensing Notice

The IEDB Analysis Resource tools (NetMHCpan, SMM, ANN, etc.) are provided under an academic/open‑source license for non‑commercial research. **For any commercial use or deployment**, you must obtain a separate commercial license:

- Visit the IEDB download page: http://tools.iedb.org/main/download/
- Follow instructions to request a commercial license.

---

## License

MIT License - See [LICENSE](LICENSE.md) for full text.