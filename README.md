# MHC-I Binding Predictor

Comprehensive tool for predicting peptide-MHC class I binding using the IEDB API, with graphical interface (GUI) and command-line interface (CLI).

## Overview

This project provides two usage modes:
1. **Modern GUI** based on PySide6 for an intuitive user experience
2. **Powerful CLI** for automation and batch processing

Predictions are made through the IEDB API using NetMHCpan methods for:
- EL binding score (Elicitation Likelihood)
- IC50 values (Binding Affinity)
- Immunogenicity calculation
- Peptide pattern analysis

## Key Features

### Graphical Interface (gui.py)
- **Direct predictions**: Input peptides and alleles to get complete scores
- **Pattern analysis**: Generation and analysis of peptide variants from patterns (e.g., `A[CD]E[FY]GH`)
- **Advanced filters**: Result filtering based on multiple thresholds
- **Export**: Save results to CSV with customizable separators
- **Tabbed interface**: Clear organization of different functionalities

### CLI Interface (mhc.py)
- **Batch predictions**: Processing peptide lists from files or direct input
- **Pattern analysis**: Variant generation and combined predictions
- **Result filtering**: Applying thresholds to existing files
- **Flexible output**: Support for different CSV separators and decimal formats

## Installation

### Prerequisites
- Python 3.7+
- Internet connection (for IEDB API)

### Dependencies
Install required dependencies:

```bash
pip install pandas numpy requests PySide6 click
```

### Project Structure
```
.
├── gui.py          # Main graphical interface
├── mhc.py          # Prediction engine and CLI
└── output/         # Results directory (created automatically)
```

## GUI Usage

Start the graphical interface:

```bash
python gui.py
```

### Available Tabs

1. **Predictions**
   - Enter peptides (one per line) or load from file
   - Specify MHC alleles (comma-separated)
   - Set peptide lengths and delay between requests
   - Click "Run Predictions" to start predictions

2. **Pattern Analysis**
   - Enter a pattern (e.g., `A[CD]E[FY]GH`)
   - Specify alleles and parameters
   - Click "Analyze Pattern" to generate variants and run predictions

3. **Filter Results**
   - Apply filters to current results
   - Available thresholds: EL Score, Percentile Rank, IC50, Immunogenicity
   - Click "Apply Filters" to update the view

### Result Features
- Sortable table with all prediction data
- Save complete or filtered results
- Export with customizable separators
- Clear view functionality

## CLI Usage

The `mhc.py` module offers a complete command-line interface:

### Available Commands

#### Predictions
```bash
python mhc.py predict --peptides "SIINFEKL,RAKFKQLL" --alleles "HLA-A*02:01,HLA-B*07:02" --delay 3.0
```

```bash
python mhc.py predict --peptides peptide_list.txt --alleles "HLA-A*02:01" --output results.csv
```

#### Pattern Analysis
```bash
python mhc.py analyze --pattern "A[CD]E[FY]GH" --alleles "HLA-A*02:01" --lengths "8,9,10" --output analysis.csv
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

## Prediction Output

Results include:
- **peptide**: Peptide sequence
- **allele**: MHC allele
- **el_score**: EL binding score (0-1, higher values indicate better binding)
- **percentile_rank**: Percentile rank (lower values indicate better binding)
- **ic50**: IC50 value in nM (lower values indicate better binding)
- **immunogenicity**: Calculated immunogenicity score

## Recommended Thresholds

For identifying strong binders:
- **EL Score**: ≥ 0
- **Percentile Rank**: ≤ 0.5
- **IC50**: ≤ 500 nM
- **Immunogenicity**: ≥ 0 (relative score)

## Technical Notes

### IEDB API
- Requests are made to `http://tools-cluster-interface.iedb.org/tools_api/mhci/`
- Methods used: `netmhcpan_el` and `netmhcpan_ba`
- Includes delay between requests to avoid API limitations

### Immunogenicity
Immunogenicity score is calculated based on:
- Specific amino acid weights
- Anchor positions for MHC allele
- Masking of key positions

### Error Handling
- The tool includes complete logging for debugging
- Network and API error handling
- Input validation

## Detailed Examples

Complete examples are available by running:
```bash
python mhc.py --help
```

## Support

For issues or questions, check:
- Internet connection
- Input format (peptides and alleles)
- Generated logs for specific errors

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
