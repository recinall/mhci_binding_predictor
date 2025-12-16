#!/usr/bin/env python3

import sys
import os
import re
import tempfile
import requests
import pandas as pd
import numpy as np
from typing import List, Dict, Set, Optional
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QTabWidget, QGroupBox, QLabel, QLineEdit, QTextEdit, QPushButton, QSpinBox,
    QDoubleSpinBox, QCheckBox, QComboBox, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QFileDialog, QMessageBox, QStatusBar, QProgressBar,
    QDialog, QDialogButtonBox, QFormLayout, QListWidget, QListWidgetItem,
    QSplitter, QFrame, QHeaderView, QAbstractScrollArea
)
from PySide6.QtCore import Qt, Signal, QThread
from PySide6.QtGui import QFont, QColor, QPalette


IMMUNOSCALE = {
    "A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110,
    "H": 0.105, "I": 0.432, "K": -0.700, "L": -0.036, "M": -0.570, "N": -0.021,
    "P": -0.036, "Q": -0.376, "R": 0.168, "S": -0.537, "T": 0.126, "V": 0.134,
    "W": 0.719, "Y": -0.012
}

IMMUNOWEIGHT = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

ALLELE_DICT = {
    "H-2-Db": "2,5,9", "H-2-Dd": "2,3,5", "H-2-Kb": "2,3,9", "H-2-Kd": "2,5,9",
    "H-2-Kk": "2,8,9", "H-2-Ld": "2,5,9", "HLA-A0101": "2,3,9", "HLA-A0201": "1,2,9",
    "HLA-A0202": "1,2,9", "HLA-A0203": "1,2,9", "HLA-A0206": "1,2,9", "HLA-A0211": "1,2,9",
    "HLA-A0301": "1,2,9", "HLA-A1101": "1,2,9", "HLA-A2301": "2,7,9", "HLA-A2402": "2,7,9",
    "HLA-A2601": "1,2,9", "HLA-A2902": "2,7,9", "HLA-A3001": "1,3,9", "HLA-A3002": "2,7,9",
    "HLA-A3101": "1,2,9", "HLA-A3201": "1,2,9", "HLA-A3301": "1,2,9", "HLA-A6801": "1,2,9",
    "HLA-A6802": "1,2,9", "HLA-A6901": "1,2,9", "HLA-B0702": "1,2,9", "HLA-B0801": "2,5,9",
    "HLA-B1501": "1,2,9", "HLA-B1502": "1,2,9", "HLA-B1801": "1,2,9", "HLA-B2705": "2,3,9",
    "HLA-B3501": "1,2,9", "HLA-B3901": "1,2,9", "HLA-B4001": "1,2,9", "HLA-B4002": "1,2,9",
    "HLA-B4402": "2,3,9", "HLA-B4403": "2,3,9", "HLA-B4501": "1,2,9", "HLA-B4601": "1,2,9",
    "HLA-B5101": "1,2,9", "HLA-B5301": "1,2,9", "HLA-B5401": "1,2,9", "HLA-B5701": "1,2,9",
    "HLA-B5801": "1,2,9"
}


def calculate_immunogenicity(peptide: str, allele: str = None) -> float:
    peptide = peptide.upper()
    peplen = len(peptide)
    cterm = peplen - 1
    score = 0.0

    if allele:
        clean_allele = allele.replace("*", "").replace(":", "")
        if clean_allele in ALLELE_DICT:
            mask_positions = [int(x) - 1 for x in ALLELE_DICT[clean_allele].split(",")]
        else:
            mask_positions = [0, 1, cterm]
    else:
        mask_positions = [0, 1, cterm]

    if peplen > 9:
        pepweight = IMMUNOWEIGHT[:5] + [0.30] * (peplen - 9) + IMMUNOWEIGHT[5:]
    else:
        pepweight = IMMUNOWEIGHT[:peplen]

    for i, aa in enumerate(peptide):
        if aa not in IMMUNOSCALE:
            return 0.0
        if i not in mask_positions and i < len(pepweight):
            score += pepweight[i] * IMMUNOSCALE[aa]

    return round(score, 5)


def tokenize_pattern(pattern: str) -> List[List[str]]:
    tokens = []
    i = 0
    while i < len(pattern):
        if pattern[i] == '[':
            close_idx = pattern.find(']', i + 1)
            if close_idx == -1:
                tokens.append([pattern[i]])
                i += 1
            else:
                options = list(pattern[i + 1:close_idx])
                tokens.append(options)
                i = close_idx + 1
        else:
            tokens.append([pattern[i]])
            i += 1
    return tokens


def generate_variants_for_length(tokens: List[List[str]], length: int) -> Set[str]:
    variants = set()
    n = len(tokens)
    if n < length:
        return variants

    for start in range(n - length + 1):
        sequence = tokens[start:start + length]

        def generate_combinations(current: str, index: int):
            if index == len(sequence):
                variants.add(current)
                return
            for option in sequence[index]:
                generate_combinations(current + option, index + 1)

        generate_combinations('', 0)

    return variants


def generate_all_variants(patterns: List[str], lengths: List[int]) -> List[str]:
    all_variants = set()
    for pattern in patterns:
        tokens = tokenize_pattern(pattern)
        for length in lengths:
            variants = generate_variants_for_length(tokens, length)
            all_variants.update(variants)
    return sorted(list(all_variants))


VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def validate_peptide(peptide: str) -> tuple:
    peptide = peptide.upper().strip()
    invalid_chars = set(peptide) - VALID_AMINO_ACIDS
    if invalid_chars:
        return False, f"Invalid characters: {', '.join(sorted(invalid_chars))}"
    if len(peptide) < 8 or len(peptide) > 15:
        return False, f"Length {len(peptide)} not in range 8-15"
    return True, peptide


class ApiWorker(QThread):
    finished = Signal(list)
    error = Signal(str)
    progress = Signal(str)

    def __init__(self, peptides: List[str], alleles: List[str], lengths: List[int], delay: float):
        super().__init__()
        self.peptides = peptides
        self.alleles = alleles
        self.lengths = lengths
        self.delay = delay
        self.api_url = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"

    def make_api_request(self, method: str, peptides: List[str], allele: str, lengths: List[int]) -> List[Dict]:
        fasta = "\n".join([f">peptide{i+1}\n{p}" for i, p in enumerate(peptides)])
        data = {
            "method": method,
            "sequence_text": fasta,
            "allele": allele,
            "length": ",".join(map(str, lengths))
        }

        try:
            response = requests.post(self.api_url, data=data, timeout=60)

            if response.status_code != 200:
                self.progress.emit(f"API error: {response.status_code}")
                return []

            text = response.text.strip()

            if "invalid character" in text.lower() or "error" in text.lower() or not text.startswith("allele\t"):
                self.progress.emit(f"API error: {text[:100]}")
                return []

            lines = text.split("\n")

            if len(lines) < 2:
                return []

            headers = lines[0].split("\t")
            results = []

            for line in lines[1:]:
                if not line.strip():
                    continue
                values = line.split("\t")
                row = {}
                for j, header in enumerate(headers):
                    if j < len(values):
                        row[header.lower()] = values[j]
                results.append(row)

            return results
        except Exception as e:
            self.progress.emit(f"Request error: {str(e)}")
            return []

    def run(self):
        try:
            valid_peptides = []
            invalid_peptides = []

            for pep in self.peptides:
                is_valid, result = validate_peptide(pep)
                if is_valid:
                    valid_peptides.append(result)
                else:
                    invalid_peptides.append(f"{pep}: {result}")

            if invalid_peptides:
                self.progress.emit(f"Skipped {len(invalid_peptides)} invalid peptides")

            if not valid_peptides:
                self.error.emit("No valid peptides to analyze")
                return

            all_results = []
            total_alleles = len(self.alleles)

            self.progress.emit(f"Analyzing {len(valid_peptides)} peptides with {total_alleles} allele(s)")

            for i, allele in enumerate(self.alleles):
                self.progress.emit(f"Processing allele {i + 1}/{total_alleles}: {allele}")

                if i > 0:
                    self.msleep(int(self.delay * 1000))

                el_results = self.make_api_request("netmhcpan_el", valid_peptides, allele, self.lengths)

                if el_results:
                    self.msleep(int(self.delay * 1000))
                    ba_results = self.make_api_request("netmhcpan_ba", valid_peptides, allele, self.lengths)

                    ba_dict = {}
                    for ba_row in ba_results:
                        pep = ba_row.get("peptide", "")
                        ba_dict[pep] = ba_row

                    for el_row in el_results:
                        combined = dict(el_row)
                        pep = el_row.get("peptide", "")
                        if pep in ba_dict:
                            ba_row = ba_dict[pep]
                            if "ic50" in ba_row:
                                combined["ic50"] = ba_row["ic50"]
                        all_results.append(combined)

            normalized = []
            for row in all_results:
                peptide = row.get("peptide") or row.get("seq") or ""
                allele_val = row.get("allele") or row.get("mhc") or ""

                el_score = None
                score_val = row.get("score") or row.get("el_score")
                if score_val is not None and score_val != "":
                    try:
                        el_score = float(score_val)
                    except (ValueError, TypeError):
                        pass

                percentile = None
                pct_val = row.get("percentile_rank") or row.get("rank")
                if pct_val is not None and pct_val != "":
                    try:
                        percentile = float(pct_val)
                    except (ValueError, TypeError):
                        pass

                ic50 = None
                ic50_val = row.get("ic50")
                if ic50_val is not None and ic50_val != "":
                    try:
                        ic50 = float(ic50_val)
                    except (ValueError, TypeError):
                        pass

                immunogenicity = calculate_immunogenicity(peptide, allele_val) if peptide else None

                normalized.append({
                    "peptide": peptide,
                    "allele": allele_val,
                    "el_score": el_score,
                    "percentile_rank": percentile,
                    "ic50": ic50,
                    "immunogenicity": immunogenicity
                })

            self.progress.emit(f"Completed: {len(normalized)} results")
            self.finished.emit(normalized)

        except Exception as e:
            self.error.emit(f"Worker error: {str(e)}")


class ResultsTable(QTableWidget):
    def __init__(self):
        super().__init__()
        self.setAlternatingRowColors(True)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.setSortingEnabled(True)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.verticalHeader().setVisible(False)
        self.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)

        self.columns = ["Peptide", "Allele", "EL Score", "Percentile Rank", "IC50 (nM)", "Immunogenicity"]
        self.setColumnCount(len(self.columns))
        self.setHorizontalHeaderLabels(self.columns)

    def load_data(self, results: List[Dict]):
        self.setSortingEnabled(False)
        self.setRowCount(len(results))

        for i, row in enumerate(results):
            peptide_val = str(row.get("peptide", ""))
            allele_val = str(row.get("allele", ""))
            self.setItem(i, 0, QTableWidgetItem(peptide_val))
            self.setItem(i, 1, QTableWidgetItem(allele_val))

            el_score = row.get("el_score")
            el_item = QTableWidgetItem(f"{el_score:.4f}" if el_score is not None else "N/A")
            el_item.setData(Qt.UserRole, el_score if el_score is not None else float('inf'))
            self.setItem(i, 2, el_item)

            percentile = row.get("percentile_rank")
            pct_item = QTableWidgetItem(f"{percentile:.2f}" if percentile is not None else "N/A")
            pct_item.setData(Qt.UserRole, percentile if percentile is not None else float('inf'))
            self.setItem(i, 3, pct_item)

            ic50 = row.get("ic50")
            ic50_item = QTableWidgetItem(f"{ic50:.2f}" if ic50 is not None else "N/A")
            ic50_item.setData(Qt.UserRole, ic50 if ic50 is not None else float('inf'))
            self.setItem(i, 4, ic50_item)

            immuno = row.get("immunogenicity")
            immuno_item = QTableWidgetItem(f"{immuno:.5f}" if immuno is not None else "N/A")
            immuno_item.setData(Qt.UserRole, immuno if immuno is not None else float('-inf'))
            self.setItem(i, 5, immuno_item)

        self.setSortingEnabled(True)

    def clear_data(self):
        self.setRowCount(0)


class PredictionTab(QWidget):
    run_prediction = Signal(list, list, list, float)

    def __init__(self):
        super().__init__()
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)

        peptides_group = QGroupBox("Peptides")
        peptides_layout = QVBoxLayout()
        self.peptides_input = QTextEdit()
        self.peptides_input.setPlaceholderText("Enter peptides (one per line)\nSIINFEKL\nRAKFKQLL\nGILGFVFTL")
        self.peptides_input.setMaximumHeight(120)
        peptides_layout.addWidget(self.peptides_input)

        load_btn = QPushButton("Load from File")
        load_btn.clicked.connect(self.load_peptides_file)
        peptides_layout.addWidget(load_btn)
        peptides_group.setLayout(peptides_layout)
        layout.addWidget(peptides_group)

        alleles_group = QGroupBox("Alleles")
        alleles_layout = QVBoxLayout()
        self.alleles_input = QLineEdit()
        self.alleles_input.setPlaceholderText("HLA-A*02:01,HLA-B*07:02")
        self.alleles_input.setText("HLA-A*02:01")
        alleles_layout.addWidget(self.alleles_input)
        alleles_group.setLayout(alleles_layout)
        layout.addWidget(alleles_group)

        params_group = QGroupBox("Parameters")
        params_layout = QGridLayout()

        params_layout.addWidget(QLabel("Lengths:"), 0, 0)
        self.lengths_input = QLineEdit("9")
        self.lengths_input.setPlaceholderText("8,9,10")
        params_layout.addWidget(self.lengths_input, 0, 1)

        params_layout.addWidget(QLabel("Delay (s):"), 1, 0)
        self.delay_spin = QDoubleSpinBox()
        self.delay_spin.setRange(0.5, 10.0)
        self.delay_spin.setValue(2.0)
        self.delay_spin.setSingleStep(0.5)
        params_layout.addWidget(self.delay_spin, 1, 1)

        params_group.setLayout(params_layout)
        layout.addWidget(params_group)

        self.run_btn = QPushButton("Run Predictions")
        self.run_btn.clicked.connect(self.on_run_clicked)
        layout.addWidget(self.run_btn)

        layout.addStretch()

    def load_peptides_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Peptides", "", "Text Files (*.txt);;All Files (*)"
        )
        if file_path:
            with open(file_path, 'r') as f:
                self.peptides_input.setPlainText(f.read())

    def on_run_clicked(self):
        peptides_text = self.peptides_input.toPlainText().strip()
        alleles_text = self.alleles_input.text().strip()
        lengths_text = self.lengths_input.text().strip()

        if not peptides_text or not alleles_text:
            QMessageBox.warning(self, "Input Error", "Please enter peptides and alleles")
            return

        peptides = [p.strip() for p in peptides_text.split('\n') if p.strip()]
        alleles = [a.strip() for a in alleles_text.split(',') if a.strip()]
        lengths = []
        for l in lengths_text.split(','):
            try:
                lengths.append(int(l.strip()))
            except ValueError:
                pass

        if not lengths:
            lengths = [9]

        self.run_prediction.emit(peptides, alleles, lengths, self.delay_spin.value())


class PatternAnalysisTab(QWidget):
    run_prediction = Signal(list, list, list, float)

    def __init__(self):
        super().__init__()
        self.current_variants = []
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)

        pattern_group = QGroupBox("Pattern Analysis")
        pattern_layout = QGridLayout()

        pattern_layout.addWidget(QLabel("Patterns:"), 0, 0)
        self.pattern_input = QLineEdit()
        self.pattern_input.setPlaceholderText("A[CD]E[FY]GH,X[LM]P[DE],Y[NQ]F")
        pattern_layout.addWidget(self.pattern_input, 0, 1)

        pattern_layout.addWidget(QLabel("Alleles:"), 1, 0)
        self.alleles_input = QLineEdit()
        self.alleles_input.setPlaceholderText("HLA-A*02:01,HLA-B*07:02")
        self.alleles_input.setText("HLA-A*02:01")
        pattern_layout.addWidget(self.alleles_input, 1, 1)

        pattern_layout.addWidget(QLabel("Lengths:"), 2, 0)
        self.lengths_input = QLineEdit("9")
        self.lengths_input.setPlaceholderText("8,9,10")
        pattern_layout.addWidget(self.lengths_input, 2, 1)

        pattern_layout.addWidget(QLabel("Delay (s):"), 3, 0)
        self.delay_spin = QDoubleSpinBox()
        self.delay_spin.setRange(0.5, 10.0)
        self.delay_spin.setValue(2.0)
        self.delay_spin.setSingleStep(0.5)
        pattern_layout.addWidget(self.delay_spin, 3, 1)

        pattern_group.setLayout(pattern_layout)
        layout.addWidget(pattern_group)

        buttons_layout = QHBoxLayout()
        self.generate_btn = QPushButton("Generate Variants")
        self.generate_btn.clicked.connect(self.generate_variants)
        buttons_layout.addWidget(self.generate_btn)

        self.analyze_btn = QPushButton("Analyze Pattern")
        self.analyze_btn.setEnabled(False)
        self.analyze_btn.clicked.connect(self.run_analysis)
        buttons_layout.addWidget(self.analyze_btn)
        layout.addLayout(buttons_layout)

        variants_group = QGroupBox("Generated Variants")
        variants_layout = QVBoxLayout()
        self.variants_display = QTextEdit()
        self.variants_display.setReadOnly(True)
        self.variants_display.setMaximumHeight(150)
        self.variants_display.setFont(QFont("Monospace", 9))
        variants_layout.addWidget(self.variants_display)
        variants_group.setLayout(variants_layout)
        layout.addWidget(variants_group)

        layout.addStretch()

    def generate_variants(self):
        patterns_text = self.pattern_input.text().strip()
        lengths_text = self.lengths_input.text().strip()

        if not patterns_text:
            QMessageBox.warning(self, "Input Error", "Please enter at least one pattern")
            return

        patterns = [p.replace(" ", "") for p in patterns_text.split(',') if p.strip()]
        lengths = []
        for l in lengths_text.split(','):
            try:
                val = int(l.strip())
                if val > 0:
                    lengths.append(val)
            except ValueError:
                pass

        if not lengths:
            lengths = [9]

        self.current_variants = generate_all_variants(patterns, lengths)
        self.variants_display.setPlainText('\n'.join(self.current_variants))
        self.analyze_btn.setEnabled(len(self.current_variants) > 0)

        QMessageBox.information(
            self, "Variants Generated",
            f"Generated {len(self.current_variants)} variants from {len(patterns)} pattern(s)"
        )

    def run_analysis(self):
        if not self.current_variants:
            QMessageBox.warning(self, "No Variants", "Please generate variants first")
            return

        alleles_text = self.alleles_input.text().strip()
        lengths_text = self.lengths_input.text().strip()

        if not alleles_text:
            QMessageBox.warning(self, "Input Error", "Please enter alleles")
            return

        alleles = [a.strip() for a in alleles_text.split(',') if a.strip()]
        lengths = []
        for l in lengths_text.split(','):
            try:
                lengths.append(int(l.strip()))
            except ValueError:
                pass

        if not lengths:
            lengths = [9]

        self.run_prediction.emit(self.current_variants, alleles, lengths, self.delay_spin.value())


class FilterTab(QWidget):
    apply_filter = Signal(dict)
    clear_filter = Signal()

    def __init__(self):
        super().__init__()
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)

        filter_group = QGroupBox("Filter Criteria")
        filter_layout = QGridLayout()

        self.peptide_check = QCheckBox("Peptide:")
        filter_layout.addWidget(self.peptide_check, 0, 0)
        self.peptide_list = QListWidget()
        self.peptide_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.peptide_list.setMaximumHeight(80)
        filter_layout.addWidget(self.peptide_list, 0, 1)

        self.allele_check = QCheckBox("Allele:")
        filter_layout.addWidget(self.allele_check, 1, 0)
        self.allele_list = QListWidget()
        self.allele_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.allele_list.setMaximumHeight(80)
        filter_layout.addWidget(self.allele_list, 1, 1)

        self.el_score_check = QCheckBox("EL Score >=")
        filter_layout.addWidget(self.el_score_check, 2, 0)
        self.el_score_spin = QDoubleSpinBox()
        self.el_score_spin.setRange(-10, 10)
        self.el_score_spin.setValue(0.0)
        self.el_score_spin.setDecimals(4)
        self.el_score_spin.setSingleStep(0.01)
        filter_layout.addWidget(self.el_score_spin, 2, 1)

        self.percentile_check = QCheckBox("Percentile <=")
        filter_layout.addWidget(self.percentile_check, 3, 0)
        self.percentile_spin = QDoubleSpinBox()
        self.percentile_spin.setRange(0, 100)
        self.percentile_spin.setValue(0.5)
        self.percentile_spin.setDecimals(2)
        filter_layout.addWidget(self.percentile_spin, 3, 1)

        self.ic50_check = QCheckBox("IC50 (nM) <=")
        filter_layout.addWidget(self.ic50_check, 4, 0)
        self.ic50_spin = QDoubleSpinBox()
        self.ic50_spin.setRange(0, 100000)
        self.ic50_spin.setValue(500)
        self.ic50_spin.setDecimals(2)
        filter_layout.addWidget(self.ic50_spin, 4, 1)

        self.immuno_check = QCheckBox("Immunogenicity >=")
        filter_layout.addWidget(self.immuno_check, 5, 0)
        self.immuno_spin = QDoubleSpinBox()
        self.immuno_spin.setRange(-10, 10)
        self.immuno_spin.setValue(0.0)
        self.immuno_spin.setDecimals(5)
        self.immuno_spin.setSingleStep(0.01)
        filter_layout.addWidget(self.immuno_spin, 5, 1)

        filter_group.setLayout(filter_layout)
        layout.addWidget(filter_group)

        buttons_layout = QHBoxLayout()
        self.apply_btn = QPushButton("Apply Filters")
        self.apply_btn.clicked.connect(self.on_apply_clicked)
        buttons_layout.addWidget(self.apply_btn)

        self.clear_btn = QPushButton("Show All Results")
        self.clear_btn.clicked.connect(self.on_clear_clicked)
        buttons_layout.addWidget(self.clear_btn)
        layout.addLayout(buttons_layout)

        export_group = QGroupBox("Export Settings")
        export_layout = QGridLayout()

        export_layout.addWidget(QLabel("CSV Separator:"), 0, 0)
        self.csv_sep_combo = QComboBox()
        self.csv_sep_combo.addItems(["Comma (,)", "Semicolon (;)", "Tab"])
        export_layout.addWidget(self.csv_sep_combo, 0, 1)

        export_layout.addWidget(QLabel("Decimal Separator:"), 1, 0)
        self.decimal_sep_combo = QComboBox()
        self.decimal_sep_combo.addItems(["Dot (.)", "Comma (,)"])
        export_layout.addWidget(self.decimal_sep_combo, 1, 1)

        export_group.setLayout(export_layout)
        layout.addWidget(export_group)

        layout.addStretch()

    def update_filter_lists(self, results: List[Dict]):
        peptides = sorted(set(r.get("peptide", "") for r in results if r.get("peptide")))
        alleles = sorted(set(r.get("allele", "") for r in results if r.get("allele")))

        self.peptide_list.clear()
        for p in peptides:
            self.peptide_list.addItem(p)

        self.allele_list.clear()
        for a in alleles:
            self.allele_list.addItem(a)

    def on_apply_clicked(self):
        filters = {}

        if self.peptide_check.isChecked():
            selected = [item.text() for item in self.peptide_list.selectedItems()]
            if selected:
                filters["peptides"] = selected

        if self.allele_check.isChecked():
            selected = [item.text() for item in self.allele_list.selectedItems()]
            if selected:
                filters["alleles"] = selected

        if self.el_score_check.isChecked():
            filters["el_score_min"] = self.el_score_spin.value()

        if self.percentile_check.isChecked():
            filters["percentile_max"] = self.percentile_spin.value()

        if self.ic50_check.isChecked():
            filters["ic50_max"] = self.ic50_spin.value()

        if self.immuno_check.isChecked():
            filters["immunogenicity_min"] = self.immuno_spin.value()

        if not filters:
            QMessageBox.warning(self, "No Filters", "Please select at least one filter criterion")
            return

        self.apply_filter.emit(filters)

    def on_clear_clicked(self):
        self.peptide_check.setChecked(False)
        self.allele_check.setChecked(False)
        self.el_score_check.setChecked(False)
        self.percentile_check.setChecked(False)
        self.ic50_check.setChecked(False)
        self.immuno_check.setChecked(False)
        self.peptide_list.clearSelection()
        self.allele_list.clearSelection()
        self.clear_filter.emit()

    def get_csv_separator(self) -> str:
        text = self.csv_sep_combo.currentText()
        if "Tab" in text:
            return "\t"
        elif "Semicolon" in text:
            return ";"
        return ","

    def get_decimal_separator(self) -> str:
        text = self.decimal_sep_combo.currentText()
        if "Comma" in text:
            return ","
        return "."


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.current_results = []
        self.filtered_results = []
        self.is_filtered = False
        self.worker = None
        self.setup_ui()

    def setup_ui(self):
        self.setWindowTitle("MHC-I Binding Predictor")
        self.setMinimumSize(1000, 700)

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)

        splitter = QSplitter(Qt.Horizontal)

        left_panel = QWidget()
        left_panel.setMaximumWidth(450)
        left_layout = QVBoxLayout(left_panel)

        self.tabs = QTabWidget()

        self.prediction_tab = PredictionTab()
        self.prediction_tab.run_prediction.connect(self.run_predictions)
        self.tabs.addTab(self.prediction_tab, "Predictions")

        self.pattern_tab = PatternAnalysisTab()
        self.pattern_tab.run_prediction.connect(self.run_predictions)
        self.tabs.addTab(self.pattern_tab, "Pattern Analysis")

        self.filter_tab = FilterTab()
        self.filter_tab.apply_filter.connect(self.apply_filters)
        self.filter_tab.clear_filter.connect(self.clear_filters)
        self.tabs.addTab(self.filter_tab, "Filter Results")

        left_layout.addWidget(self.tabs)
        splitter.addWidget(left_panel)

        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)

        header_layout = QHBoxLayout()
        results_label = QLabel("Results")
        results_label.setFont(QFont("Arial", 14, QFont.Bold))
        header_layout.addWidget(results_label)
        header_layout.addStretch()

        export_btn = QPushButton("Export CSV")
        export_btn.clicked.connect(self.export_results)
        header_layout.addWidget(export_btn)

        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_results)
        header_layout.addWidget(clear_btn)

        right_layout.addLayout(header_layout)

        self.results_table = ResultsTable()
        right_layout.addWidget(self.results_table)

        splitter.addWidget(right_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        main_layout.addWidget(splitter)

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setMaximumWidth(200)
        self.status_bar.addPermanentWidget(self.progress_bar)

        self.results_count_label = QLabel("0 results")
        self.status_bar.addPermanentWidget(self.results_count_label)

        self.status_bar.showMessage("Ready")

    def run_predictions(self, peptides: List[str], alleles: List[str], lengths: List[int], delay: float):
        if self.worker and self.worker.isRunning():
            QMessageBox.warning(self, "Busy", "Another operation is in progress")
            return

        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.set_controls_enabled(False)

        self.worker = ApiWorker(peptides, alleles, lengths, delay)
        self.worker.finished.connect(self.on_predictions_finished)
        self.worker.error.connect(self.on_predictions_error)
        self.worker.progress.connect(self.on_progress_update)
        self.worker.start()

    def on_predictions_finished(self, results: List[Dict]):
        self.progress_bar.setVisible(False)
        self.set_controls_enabled(True)

        self.current_results = results
        self.filtered_results = []
        self.is_filtered = False

        self.results_table.load_data(results)
        self.filter_tab.update_filter_lists(results)
        self.results_count_label.setText(f"{len(results)} results")
        self.status_bar.showMessage(f"Predictions completed: {len(results)} results")

    def on_predictions_error(self, error_msg: str):
        self.progress_bar.setVisible(False)
        self.set_controls_enabled(True)
        QMessageBox.critical(self, "Error", f"Prediction failed: {error_msg}")
        self.status_bar.showMessage("Error occurred")

    def on_progress_update(self, message: str):
        self.status_bar.showMessage(message)

    def set_controls_enabled(self, enabled: bool):
        self.prediction_tab.run_btn.setEnabled(enabled)
        self.pattern_tab.generate_btn.setEnabled(enabled)
        self.pattern_tab.analyze_btn.setEnabled(enabled and len(self.pattern_tab.current_variants) > 0)
        self.filter_tab.apply_btn.setEnabled(enabled)

    def apply_filters(self, filters: Dict):
        if not self.current_results:
            QMessageBox.warning(self, "No Data", "No results to filter")
            return

        filtered = list(self.current_results)

        if "peptides" in filters:
            filtered = [r for r in filtered if r.get("peptide") in filters["peptides"]]

        if "alleles" in filters:
            filtered = [r for r in filtered if r.get("allele") in filters["alleles"]]

        if "el_score_min" in filters:
            threshold = filters["el_score_min"]
            filtered = [r for r in filtered if r.get("el_score") is not None and r["el_score"] >= threshold]

        if "percentile_max" in filters:
            threshold = filters["percentile_max"]
            filtered = [r for r in filtered if r.get("percentile_rank") is not None and r["percentile_rank"] <= threshold]

        if "ic50_max" in filters:
            threshold = filters["ic50_max"]
            filtered = [r for r in filtered if r.get("ic50") is not None and r["ic50"] <= threshold]

        if "immunogenicity_min" in filters:
            threshold = filters["immunogenicity_min"]
            filtered = [r for r in filtered if r.get("immunogenicity") is not None and r["immunogenicity"] >= threshold]

        self.filtered_results = filtered
        self.is_filtered = True
        self.results_table.load_data(filtered)
        self.results_count_label.setText(f"{len(filtered)} results (filtered)")
        self.status_bar.showMessage(f"Filtered to {len(filtered)} results")

    def clear_filters(self):
        if not self.current_results:
            return

        self.is_filtered = False
        self.filtered_results = []
        self.results_table.load_data(self.current_results)
        self.results_count_label.setText(f"{len(self.current_results)} results")
        self.status_bar.showMessage("Showing all results")

    def export_results(self):
        results = self.filtered_results if self.is_filtered else self.current_results

        if not results:
            QMessageBox.warning(self, "No Data", "No results to export")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results", f"mhc_binding_results.csv", "CSV Files (*.csv)"
        )

        if not file_path:
            return

        csv_sep = self.filter_tab.get_csv_separator()
        decimal_sep = self.filter_tab.get_decimal_separator()

        def format_number(val, decimals):
            if val is None:
                return ""
            formatted = f"{val:.{decimals}f}"
            if decimal_sep == ",":
                formatted = formatted.replace(".", ",")
            return formatted

        headers = ["Peptide", "Allele", "EL Score", "Percentile Rank", "IC50 (nM)", "Immunogenicity"]
        lines = [csv_sep.join(headers)]

        for r in results:
            row = [
                r.get("peptide", ""),
                r.get("allele", ""),
                format_number(r.get("el_score"), 4),
                format_number(r.get("percentile_rank"), 2),
                format_number(r.get("ic50"), 2),
                format_number(r.get("immunogenicity"), 5)
            ]
            lines.append(csv_sep.join(row))

        with open(file_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))

        result_type = "filtered" if self.is_filtered else "all"
        self.status_bar.showMessage(f"Exported {len(results)} {result_type} results to {file_path}")

    def clear_results(self):
        self.current_results = []
        self.filtered_results = []
        self.is_filtered = False
        self.results_table.clear_data()
        self.filter_tab.peptide_list.clear()
        self.filter_tab.allele_list.clear()
        self.results_count_label.setText("0 results")
        self.status_bar.showMessage("Results cleared")


def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')

    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(245, 245, 245))
    palette.setColor(QPalette.WindowText, QColor(44, 62, 80))
    palette.setColor(QPalette.Base, QColor(255, 255, 255))
    palette.setColor(QPalette.AlternateBase, QColor(248, 249, 250))
    palette.setColor(QPalette.ToolTipBase, QColor(255, 255, 255))
    palette.setColor(QPalette.ToolTipText, QColor(44, 62, 80))
    palette.setColor(QPalette.Text, QColor(44, 62, 80))
    palette.setColor(QPalette.Button, QColor(255, 255, 255))
    palette.setColor(QPalette.ButtonText, QColor(44, 62, 80))
    palette.setColor(QPalette.Highlight, QColor(102, 126, 234))
    palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))
    app.setPalette(palette)

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()
