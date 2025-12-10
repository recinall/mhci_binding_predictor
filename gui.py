#!/usr/bin/env python3
"""
GUI for MHC-I Binding Predictor using PySide6
Modern interface for peptide-MHC binding prediction using IEDB API
"""

import sys
import os
import threading
from pathlib import Path
from PySide6.QtWidgets import *
from PySide6.QtCore import *
from PySide6.QtGui import *
import pandas as pd

from mhc import IEDBBindingPredictor

class WorkerThread(QThread):
    finished = Signal(pd.DataFrame)
    error = Signal(str)
    progress = Signal(str)

    def __init__(self, predictor, operation, **kwargs):
        super().__init__()
        self.predictor = predictor
        self.operation = operation
        self.kwargs = kwargs

    def run(self):
        try:
            if self.operation == "predict":
                self.progress.emit("Starting predictions...")
                result = self.predictor.predict_comprehensive(
                    self.kwargs['peptides'], 
                    self.kwargs['alleles'], 
                    self.kwargs['lengths'],
                    delay=self.kwargs.get('delay', 2.0)
                )
                self.finished.emit(result)
            elif self.operation == "analyze":
                self.progress.emit("Generating variants...")
                variants = self.predictor.generate_variants(self.kwargs['pattern'])
                self.progress.emit("Running predictions...")
                result = self.predictor.predict_comprehensive(
                    variants,
                    self.kwargs['alleles'],
                    self.kwargs['lengths'],
                    delay=self.kwargs.get('delay', 2.0)
                )
                self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))

class ResultsTable(QTableWidget):
    def __init__(self):
        super().__init__()
        self.setAlternatingRowColors(True)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setSortingEnabled(True)
        self.horizontalHeader().setStretchLastSection(True)

    def load_data(self, df):
        if df.empty:
            return
        
        self.setRowCount(len(df))
        self.setColumnCount(len(df.columns))
        self.setHorizontalHeaderLabels(df.columns.tolist())
        
        for i, row in df.iterrows():
            for j, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                self.setItem(i, j, item)
        
        self.resizeColumnsToContents()

class PredictionWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.predictor = parent.predictor
        self.parent = parent
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        peptides_group = QGroupBox("Peptides")
        peptides_layout = QVBoxLayout()
        
        self.peptides_text = QTextEdit()
        self.peptides_text.setPlaceholderText("Enter peptides (one per line) or click 'Load from file'")
        self.peptides_text.setMaximumHeight(100)
        peptides_layout.addWidget(self.peptides_text)
        
        self.load_peptides_btn = QPushButton("Load from file")
        self.load_peptides_btn.clicked.connect(self.load_peptides_file)
        peptides_layout.addWidget(self.load_peptides_btn)
        
        peptides_group.setLayout(peptides_layout)
        layout.addWidget(peptides_group)

        alleles_group = QGroupBox("Alleles")
        alleles_layout = QVBoxLayout()
        
        self.alleles_text = QLineEdit()
        self.alleles_text.setPlaceholderText("HLA-A*02:01,HLA-B*07:02")
        alleles_layout.addWidget(self.alleles_text)
        
        alleles_group.setLayout(alleles_layout)
        layout.addWidget(alleles_group)

        params_group = QGroupBox("Parameters")
        params_layout = QGridLayout()
        
        params_layout.addWidget(QLabel("Lengths:"), 0, 0)
        self.lengths_text = QLineEdit("9")
        params_layout.addWidget(self.lengths_text, 0, 1)
        
        params_layout.addWidget(QLabel("Delay (s):"), 1, 0)
        self.delay_spin = QDoubleSpinBox()
        self.delay_spin.setRange(0.1, 10.0)
        self.delay_spin.setValue(2.0)
        params_layout.addWidget(self.delay_spin, 1, 1)
        
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)

        self.predict_btn = QPushButton("Run Predictions")
        self.predict_btn.clicked.connect(self.run_predictions)
        layout.addWidget(self.predict_btn)

        self.setLayout(layout)

    def load_peptides_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Peptides", "", "Text Files (*.txt)")
        if file_path:
            with open(file_path, 'r') as f:
                peptides = f.read()
            self.peptides_text.setPlainText(peptides)

    def run_predictions(self):
        peptides_text = self.peptides_text.toPlainText().strip()
        if not peptides_text or not self.alleles_text.text().strip():
            QMessageBox.warning(self, "Input Error", "Please enter peptides and alleles")
            return

        peptides = [p.strip() for p in peptides_text.split('\n') if p.strip()]
        alleles = [a.strip() for a in self.alleles_text.text().split(',') if a.strip()]
        lengths = [int(l.strip()) for l in self.lengths_text.text().split(',') if l.strip()]

        self.parent.run_worker_thread("predict", 
                                    peptides=peptides, 
                                    alleles=alleles, 
                                    lengths=lengths,
                                    delay=self.delay_spin.value())

class AnalysisWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.predictor = parent.predictor
        self.parent = parent
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        pattern_group = QGroupBox("Pattern Analysis")
        pattern_layout = QGridLayout()
        
        pattern_layout.addWidget(QLabel("Pattern:"), 0, 0)
        self.pattern_text = QLineEdit()
        self.pattern_text.setPlaceholderText("A[CD]E[FY]GH")
        pattern_layout.addWidget(self.pattern_text, 0, 1)
        
        pattern_layout.addWidget(QLabel("Alleles:"), 1, 0)
        self.alleles_text = QLineEdit()
        self.alleles_text.setPlaceholderText("HLA-A*02:01,HLA-B*07:02")
        pattern_layout.addWidget(self.alleles_text, 1, 1)
        
        pattern_layout.addWidget(QLabel("Lengths:"), 2, 0)
        self.lengths_text = QLineEdit("9")
        pattern_layout.addWidget(self.lengths_text, 2, 1)
        
        pattern_layout.addWidget(QLabel("Delay (s):"), 3, 0)
        self.delay_spin = QDoubleSpinBox()
        self.delay_spin.setRange(0.1, 10.0)
        self.delay_spin.setValue(2.0)
        pattern_layout.addWidget(self.delay_spin, 3, 1)
        
        pattern_group.setLayout(pattern_layout)
        layout.addWidget(pattern_group)

        self.analyze_btn = QPushButton("Analyze Pattern")
        self.analyze_btn.clicked.connect(self.run_analysis)
        layout.addWidget(self.analyze_btn)

        variants_group = QGroupBox("Generated Variants")
        variants_layout = QVBoxLayout()
        
        self.variants_text = QTextEdit()
        self.variants_text.setMaximumHeight(150)
        self.variants_text.setReadOnly(True)
        variants_layout.addWidget(self.variants_text)
        
        variants_group.setLayout(variants_layout)
        layout.addWidget(variants_group)

        self.setLayout(layout)

    def run_analysis(self):
        pattern = self.pattern_text.text().replace(" ", "")
        alleles_text = self.alleles_text.text().strip()
        
        if not pattern or not alleles_text:
            QMessageBox.warning(self, "Input Error", "Please enter pattern and alleles")
            return

        variants = self.predictor.generate_variants(pattern)
        self.variants_text.setPlainText('\n'.join(variants))

        alleles = [a.strip() for a in alleles_text.split(',') if a.strip()]
        lengths = [int(l.strip()) for l in self.lengths_text.text().split(',') if l.strip()]

        self.parent.run_worker_thread("analyze", 
                                    pattern=pattern,
                                    alleles=alleles, 
                                    lengths=lengths,
                                    delay=self.delay_spin.value())

class FilterWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        filter_group = QGroupBox("Filter Parameters")
        filter_layout = QGridLayout()

        self.el_score_check = QCheckBox("EL Score >=")
        filter_layout.addWidget(self.el_score_check, 0, 0)
        self.el_score_spin = QDoubleSpinBox()
        self.el_score_spin.setRange(-10, 10)
        self.el_score_spin.setValue(0.5)
        filter_layout.addWidget(self.el_score_spin, 0, 1)

        self.percentile_check = QCheckBox("Percentile <=")
        filter_layout.addWidget(self.percentile_check, 1, 0)
        self.percentile_spin = QDoubleSpinBox()
        self.percentile_spin.setRange(0, 100)
        self.percentile_spin.setValue(2.0)
        filter_layout.addWidget(self.percentile_spin, 1, 1)

        self.ic50_check = QCheckBox("IC50 <=")
        filter_layout.addWidget(self.ic50_check, 2, 0)
        self.ic50_spin = QDoubleSpinBox()
        self.ic50_spin.setRange(0, 10000)
        self.ic50_spin.setValue(500)
        filter_layout.addWidget(self.ic50_spin, 2, 1)

        self.immuno_check = QCheckBox("Immunogenicity >=")
        filter_layout.addWidget(self.immuno_check, 3, 0)
        self.immuno_spin = QDoubleSpinBox()
        self.immuno_spin.setRange(-10, 10)
        self.immuno_spin.setValue(0.5)
        filter_layout.addWidget(self.immuno_spin, 3, 1)

        filter_group.setLayout(filter_layout)
        layout.addWidget(filter_group)

        self.filter_btn = QPushButton("Apply Filters")
        self.filter_btn.clicked.connect(self.apply_filters)
        layout.addWidget(self.filter_btn)

        self.setLayout(layout)

    def apply_filters(self):
        if self.parent.current_results is None or self.parent.current_results.empty:
            QMessageBox.warning(self, "No Data", "No results to filter")
            return

        kwargs = {}
        if self.el_score_check.isChecked():
            kwargs['el_score_threshold'] = self.el_score_spin.value()
        if self.percentile_check.isChecked():
            kwargs['percentile_threshold'] = self.percentile_spin.value()
        if self.ic50_check.isChecked():
            kwargs['ic50_threshold'] = self.ic50_spin.value()
        if self.immuno_check.isChecked():
            kwargs['immunogenicity_threshold'] = self.immuno_spin.value()

        if not kwargs:
            QMessageBox.warning(self, "No Filters", "Please select at least one filter")
            return

        filtered_df = self.parent.predictor.filter_binders(self.parent.current_results, **kwargs)
        self.parent.filtered_results = filtered_df
        self.parent.is_filtered_view = True
        self.parent.results_table.load_data(filtered_df)
        self.parent.status_bar.showMessage(f"Filtered to {len(filtered_df)} results")

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.predictor = IEDBBindingPredictor()
        self.current_results = None
        self.filtered_results = None
        self.is_filtered_view = False
        self.worker_thread = None
        self.setup_ui()
        
    def setup_ui(self):
        self.setWindowTitle("MHC-I Binding Predictor")
        self.setMinimumSize(800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        layout = QHBoxLayout()
        
        left_panel = QWidget()
        left_panel.setMaximumWidth(400)
        left_layout = QVBoxLayout()

        self.tab_widget = QTabWidget()
        
        self.prediction_widget = PredictionWidget(self)
        self.tab_widget.addTab(self.prediction_widget, "Predictions")
        
        self.analysis_widget = AnalysisWidget(self)
        self.tab_widget.addTab(self.analysis_widget, "Pattern Analysis")
        
        self.filter_widget = FilterWidget(self)
        self.tab_widget.addTab(self.filter_widget, "Filter Results")
        
        left_layout.addWidget(self.tab_widget)
        left_panel.setLayout(left_layout)
        layout.addWidget(left_panel)

        right_panel = QWidget()
        right_layout = QVBoxLayout()
        
        results_label = QLabel("Results")
        results_label.setFont(QFont("Arial", 12, QFont.Bold))
        right_layout.addWidget(results_label)
        
        self.results_table = ResultsTable()
        right_layout.addWidget(self.results_table)
        
        buttons_layout = QHBoxLayout()
        
        self.save_btn = QPushButton("Save Results")
        self.save_btn.clicked.connect(self.save_results)
        buttons_layout.addWidget(self.save_btn)
        
        self.clear_btn = QPushButton("Clear Results")
        self.clear_btn.clicked.connect(self.clear_results)
        buttons_layout.addWidget(self.clear_btn)

        self.exSet_btn = QPushButton("Export Settings")
        self.exSet_btn.clicked.connect(self.show_export_settings)
        buttons_layout.addWidget(self.exSet_btn)

        right_layout.addLayout(buttons_layout)
        right_panel.setLayout(right_layout)
        layout.addWidget(right_panel)

        central_widget.setLayout(layout)

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.status_bar.addPermanentWidget(self.progress_bar)

    def run_worker_thread(self, operation, **kwargs):
        if self.worker_thread and self.worker_thread.isRunning():
            QMessageBox.warning(self, "Busy", "Another operation is running")
            return

        self.worker_thread = WorkerThread(self.predictor, operation, **kwargs)
        self.worker_thread.finished.connect(self.on_worker_finished)
        self.worker_thread.error.connect(self.on_worker_error)
        self.worker_thread.progress.connect(self.on_worker_progress)
        
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.worker_thread.start()

    def on_worker_finished(self, result):
        self.progress_bar.setVisible(False)
        self.current_results = result
        self.is_filtered_view = False
        self.filtered_results = None
        self.results_table.load_data(result)
        self.status_bar.showMessage(f"Completed. {len(result)} results")

    def on_worker_error(self, error_msg):
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Error", f"Operation failed: {error_msg}")
        self.status_bar.showMessage("Error occurred")

    def on_worker_progress(self, message):
        self.status_bar.showMessage(message)

    def show_export_settings(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Export Settings")
        dialog.setFixedSize(300, 150)
        
        layout = QFormLayout()
        
        csv_sep_combo = QComboBox()
        csv_sep_combo.addItems([",", ";", "\t"])
        csv_sep_combo.setCurrentText(self.predictor.csv_separator)
        layout.addRow("CSV Separator:", csv_sep_combo)
        
        decimal_sep_combo = QComboBox()
        decimal_sep_combo.addItems([".", ","])
        decimal_sep_combo.setCurrentText(self.predictor.decimal_separator)
        layout.addRow("Decimal Separator:", decimal_sep_combo)
        
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)
        
        dialog.setLayout(layout)
        
        if dialog.exec() == QDialog.Accepted:
            self.predictor.csv_separator = csv_sep_combo.currentText()
            self.predictor.decimal_separator = decimal_sep_combo.currentText()
            self.status_bar.showMessage(f"Export settings updated: sep='{self.predictor.csv_separator}', decimal='{self.predictor.decimal_separator}'")

    def save_results(self):
        if self.current_results is None or self.current_results.empty:
            QMessageBox.warning(self, "No Data", "No results to save")
            return

        file_path, _ = QFileDialog.getSaveFileName(self, "Save Results", "", "CSV Files (*.csv)")
        if file_path:
            data_to_save = self.filtered_results if self.is_filtered_view and self.filtered_results is not None else self.current_results
            self.predictor.save_to_csv(data_to_save, file_path)
            view_type = "filtered" if self.is_filtered_view and self.filtered_results is not None else "all"
            self.status_bar.showMessage(f"Results ({view_type}) saved to {file_path}")

    def clear_results(self):
        self.current_results = None
        self.filtered_results = None
        self.is_filtered_view = False
        self.results_table.clear()
        self.results_table.setRowCount(0)
        self.results_table.setColumnCount(0)
        self.status_bar.showMessage("Results cleared")

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
