# MHC-I Binding Predictor

Comprehensive tool for predicting peptide-MHC class I binding using the IEDB API, with graphical interfaces (HTML/Python GUI) and command-line interface (CLI).

## Overview

This project provides three usage modes:
1. **Web GUI** (gui.html) - Standalone HTML interface that runs in any browser
2. **Desktop GUI** (gui.py) - Native Python application based on PySide6
3. **CLI** (mhc.py) - Command-line interface for automation and batch processing

Predictions are made through the IEDB API using NetMHCpan methods for:
- EL binding score (Elicitation Likelihood)
- IC50 values (Binding Affinity)
- Immunogenicity calculation
- Peptide pattern analysis with variant generation

## Key Features

### Web Interface (gui.html)
- **No installation required**: Works directly in any modern browser
- **Direct predictions**: Input peptides and alleles to get complete scores
- **Multi-pattern analysis**: Support for multiple patterns separated by commas
- **Advanced filters**: Filter by peptide, allele, EL Score, Percentile, IC50, Immunogenicity
- **Export settings**: Customizable CSV and decimal separators
- **Responsive design**: Works on desktop and mobile devices

### Desktop GUI (gui.py)
- **Native application**: Fast and responsive PySide6 interface
- **Direct predictions**: Input peptides and alleles to get complete scores
- **Pattern analysis**: Generation and analysis of peptide variants from patterns
- **Advanced filters**: Result filtering with peptide/allele selection lists
- **Sortable results**: Click column headers to sort data
- **Export**: Save results to CSV with customizable separators

### CLI Interface (mhc.py)
- **Batch predictions**: Processing peptide lists from files or direct input
- **Pattern analysis**: Variant generation and combined predictions
- **Result filtering**: Applying thresholds to existing files
- **Flexible output**: Support for different CSV separators and decimal formats

## Installation

### Prerequisites
- Python 3.8+
- Internet connection (for IEDB API)

### Dependencies
Install required dependencies:

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install pandas numpy requests PySide6 click
```

### Project Structure
```
.
├── gui.html        # Web interface (standalone)
├── gui.py          # Desktop GUI application
├── mhc.py          # Prediction engine and CLI
├── requirements.txt
└── output/         # Results directory (created automatically)
```

## Web GUI Usage

Simply open `gui.html` in a web browser:

```bash
firefox gui.html
# or
google-chrome gui.html
```

**Note**: For full functionality, serve the file via a local web server:

```bash
python3 -m http.server 8000
# Then open: http://localhost:8000/gui.html
```

### Tabs

1. **Predictions**
   - Enter peptides (one per line)
   - Specify MHC alleles (comma-separated)
   - Set peptide lengths and delay between requests
   - Click "Run Predictions"

2. **Pattern Analysis**
   - Enter patterns (comma-separated, e.g., `A[CD]EFGHIK,X[LM]PQRST`)
   - Specify alleles and lengths
   - Click "Generate Variants" then "Analyze Pattern"

3. **Filter Results**
   - Select peptides/alleles from dropdown lists
   - Set numeric thresholds (EL Score, Percentile, IC50, Immunogenicity)
   - Configure export settings (CSV/decimal separators)

## Desktop GUI Usage

Start the graphical interface:

```bash
python gui.py
```

### Tabs

1. **Predictions**
   - Enter peptides (one per line) or load from file
   - Specify MHC alleles (comma-separated)
   - Set peptide lengths and delay between requests
   - Click "Run Predictions"

2. **Pattern Analysis**
   - Enter patterns (comma-separated)
   - Specify alleles and parameters
   - Click "Generate Variants" then "Analyze Pattern"

3. **Filter Results**
   - Select peptides/alleles from lists (multi-select supported)
   - Set numeric thresholds
   - Configure export settings

### Result Features
- Sortable table with all prediction data
- Save complete or filtered results
- Export with customizable separators

## CLI Usage

The `mhc.py` module offers a complete command-line interface:

### Available Commands

#### Predictions
```bash
# Direct peptide input
python mhc.py predict --peptides "SIINFEKL,RAKFKQLL" --alleles "HLA-A*02:01" --lengths "8" --delay 2.0

# From file
python mhc.py predict --peptides peptide_list.txt --alleles "HLA-A*02:01,HLA-B*07:02" --output results.csv
```

#### Pattern Analysis
```bash
python mhc.py analyze --pattern "A[CD]E[FY]GHI" --alleles "HLA-A*02:01" --lengths "8,9" --output analysis.csv
```

#### Filtering
```bash
python mhc.py filter --input results.csv --percentile 2.0 --ic50 500 --immunogenicity 0.5
```

#### Variant Generation
```bash
python mhc.py variants "A[CD]E[FY]GH" --output variants.txt
```

### Global Options
- `--output-dir`: Output directory (default: `./output`)
- `--csv-sep`: CSV column separator (default: `,`)
- `--decimal-sep`: Decimal separator (default: `.`)

## Pattern Syntax

Patterns use bracket notation to specify variable positions:
- `A[CD]E` generates: `ACE`, `ADE`
- `[AB][XY]` generates: `AX`, `AY`, `BX`, `BY`

Multiple patterns can be specified (comma-separated in GUI, single pattern in CLI).

**Important**: The `lengths` parameter determines which subsequences are extracted from the pattern. For a pattern of 10 characters with `lengths=9`, two 9-mer variants are generated per combination.

## Prediction Output

Results include:
- **peptide**: Peptide sequence
- **allele**: MHC allele
- **el_score**: EL binding score (0-1, higher = better binding)
- **percentile_rank**: Percentile rank (lower = better binding)
- **ic50**: IC50 value in nM (lower = better binding)
- **immunogenicity**: Calculated immunogenicity score

## Recommended Thresholds

For identifying strong binders:

| Metric | Strong Binder | Weak Binder |
|--------|---------------|-------------|
| Percentile Rank | ≤ 0.5 | ≤ 2.0 |
| IC50 | ≤ 50 nM | ≤ 500 nM |
| EL Score | ≥ 0.5 | ≥ 0.1 |

## Technical Notes

### IEDB API
- Endpoint: `https://tools-cluster-interface.iedb.org/tools_api/mhci/`
- Methods: `netmhcpan_el` (elution) and `netmhcpan_ba` (binding affinity)
- Configurable delay between requests (default: 2 seconds)

### Immunogenicity Calculation
Based on:
- Amino acid immunogenicity weights
- Allele-specific anchor positions
- Position-dependent masking

### Length Parameter
The `length` parameter must match the peptide length for predictions to work:
- 8-mer peptides require `length=8`
- 9-mer peptides require `length=9`
- For mixed lengths, specify all: `length=8,9,10`

## Troubleshooting

### No results returned
1. Check that `lengths` matches your peptide lengths
2. Verify allele format (e.g., `HLA-A*02:01`)
3. Check internet connection

### CORS errors (Web GUI)
Serve the file via a local web server instead of opening directly:
```bash
python3 -m http.server 8000
```

### API timeout
Increase the delay between requests to avoid rate limiting.

## IEDB Data Attribution & License

This toolkit utilizes data from the Immune Epitope Database (IEDB). All IEDB data are released under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

1. **Credit IEDB**: Include attribution in publications:
   *Data courtesy of the Immune Epitope Database (IEDB; https://www.iedb.org) under CC BY 4.0.*
2. **License**: https://creativecommons.org/licenses/by/4.0/

## API Rate Limiting

- **Throttle Requests**: Default 2-second delay between API calls
- **Handle 429 Responses**: Implement back-off on rate limit errors
- **Parallelism**: Requests are made sequentially per allele

## Commercial Licensing Notice

The IEDB Analysis Resource tools are provided under an academic/open-source license for non-commercial research. **For commercial use**, obtain a separate license:

- Visit: http://tools.iedb.org/main/download/

## License

MIT License - See [LICENSE](LICENSE.md) for full text.
