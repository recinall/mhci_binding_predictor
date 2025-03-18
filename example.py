#!/usr/bin/env python3
"""
Example of using the peptide_analysis library.
"""

import os
import argparse
from peptide_analysis import (
    generate_and_analyze_peptides,
    run_complete_analysis,
    load_results_from_csv,
    filter_results_by_percentile,
    plot_score_distribution,
    plot_percentile_distribution,
    plot_allele_comparison
)

def main():
    parser = argparse.ArgumentParser(description="Analysis of 9-mer peptides")
    
    # General arguments
    parser.add_argument("--mode", choices=["generate", "analyze", "visualize"], default="analyze",
                        help="Execution mode")
    parser.add_argument("--output-dir", default="output",
                        help="Output directory")
    
    # Arguments for generation
    parser.add_argument("--num-peptides", type=int, default=1000,
                        help="Number of peptides to generate")
    
    # Arguments for analysis
    parser.add_argument("--input-csv", 
                        help="Input CSV file with peptides")
    parser.add_argument("--alleles", default="HLA-A*01:01,HLA-A*02:01",
                        help="Comma-separated list of HLA alleles")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Batch size for requests")
    parser.add_argument("--percentile-threshold", type=float, default=10.0,
                        help="Percentile rank threshold for filtering")
    parser.add_argument("--percentile-operator", default="<=", choices=["<", "<=", ">", ">=", "=="],
                        help="Comparison operator for percentile rank filtering")
    
    # Arguments for visualization
    parser.add_argument("--results-csv", 
                        help="CSV file with analysis results")
    
    args = parser.parse_args()
    
    # Execute the requested operation
    if args.mode == "generate":
        # Generate and analyze peptides
        peptides, results = generate_and_analyze_peptides(
            args.num_peptides, args.alleles, args.batch_size, args.output_dir
        )
        print(f"Generated and analyzed {len(peptides)} peptides")
        
    elif args.mode == "analyze":
        # Perform complete analysis
        report = run_complete_analysis(
            args.input_csv, args.num_peptides, args.alleles, 
            args.batch_size, args.output_dir, args.percentile_threshold,
            args.percentile_operator
        )
        print("\nAnalysis report:")
        for key, value in report.items():
            print(f"- {key}: {value}")
        
    elif args.mode == "visualize":
        # Visualize the results
        if not args.results_csv:
            print("Error: you must specify --results-csv for visualize mode")
            return
        
        results = load_results_from_csv(args.results_csv)
        if not results:
            print(f"No results found in {args.results_csv}")
            return
            
        # Filter results if a threshold was specified
        if args.percentile_threshold is not None:
            results = filter_results_by_percentile(
                results, args.percentile_threshold, args.percentile_operator
            )
        
        plots_dir = os.path.join(args.output_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)
        
        plot_score_distribution(results, os.path.join(plots_dir, "score_distribution.png"))
        plot_percentile_distribution(results, os.path.join(plots_dir, "percentile_distribution.png"))
        plot_allele_comparison(results, os.path.join(plots_dir, "allele_comparison.png"))
        
        print(f"Graphs saved to {plots_dir}")

if __name__ == "__main__":
    main()
