from setuptools import setup, find_packages

setup(
    name="peptide_analysis",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "requests",
        "pandas",
        "matplotlib",
        "seaborn",
        "tqdm",
        "numpy",
        "openpyxl",  # Per il supporto Excel
    ],
    author="Peptide Analysis Team",
    author_email="info@peptideanalysis.org",
    description="Una libreria per la generazione, analisi e visualizzazione di peptidi 9-mer",
    keywords="peptide, bioinformatics, MHC, HLA, immunology",
    python_requires=">=3.6",
)
