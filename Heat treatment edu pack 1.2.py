import sys
import os
import csv
import json
import re
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from scipy.optimize import root, root_scalar

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QGroupBox, QComboBox, QLineEdit, QTextEdit,
    QTableWidget, QTableWidgetItem, QHeaderView, QFileDialog, QMessageBox,
    QTabWidget, QDoubleSpinBox, QCheckBox, QDialog, QDialogButtonBox,
    QSplitter, QFrame, QSizePolicy, QScrollArea, QFormLayout, QGridLayout,
    QStyle
)
from PyQt5.QtCore import Qt, QFile, QFileInfo
from PyQt5.QtGui import QFont, QPalette, QColor, QIcon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.backend_bases import MouseEvent

# =============================================================================
# Constants
# =============================================================================
R = 8.314459
K = 273.15

# =============================================================================
# 1. TEMPERDATA – Steel Tempering Analysis (zwiększona czcionka, ikony)
# =============================================================================
class TemperingCurveApp(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Temperdata v 1.2 - Steel Tempering Analysis")
        self.setGeometry(100, 100, 1200, 800)
        # Zwiększenie czcionki
        font = QFont("Segoe UI", 10)
        self.setFont(font)
        
        self.data = None
        self.current_steel_type = None
        self.current_steel_data = None
        self.annotation = None
        self.plots = []
        self.plot_data = []
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        self.tab1 = QWidget()
        self.tabs.addTab(self.tab1, "Plot")
        self.tab2 = QWidget()
        self.tabs.addTab(self.tab2, "Steel Information")
        self.tab3 = QWidget()
        self.tabs.addTab(self.tab3, "Data Table")
        self.tab4 = QWidget()
        self.tabs.addTab(self.tab4, "Description & Manual")
        self.init_plot_tab()
        self.init_info_tab()
        self.init_table_tab()
        self.init_manual_tab()
        self.initialize_plot()
        self.try_auto_load_data()

    def init_plot_tab(self):
        layout = QVBoxLayout(self.tab1)
        controls_layout = QHBoxLayout()
        self.steel_label = QLabel("Select steel grade:")
        self.steel_combo = QComboBox()
        self.steel_combo.currentTextChanged.connect(self.update_all)
        self.load_button = QPushButton(" Load CSV Data")
        self.load_button.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.load_button.clicked.connect(self.load_data)
        self.save_plot_button = QPushButton(" Save Plot")
        self.save_plot_button.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_plot_button.clicked.connect(self.save_plot)
        self.save_plot_button.setEnabled(False)
        controls_layout.addWidget(self.steel_label)
        controls_layout.addWidget(self.steel_combo)
        controls_layout.addWidget(self.load_button)
        controls_layout.addWidget(self.save_plot_button)
        self.figure, self.ax = plt.subplots(figsize=(10, 6))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        layout.addLayout(controls_layout)
        layout.addWidget(self.canvas)

    def init_info_tab(self):
        layout = QVBoxLayout(self.tab2)
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setFont(QFont("Segoe UI", 10))
        layout.addWidget(self.info_text)

    def init_table_tab(self):
        layout = QVBoxLayout(self.tab3)
        button_layout = QHBoxLayout()
        self.save_table_button = QPushButton(" Save Table to CSV")
        self.save_table_button.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_table_button.clicked.connect(self.save_table)
        self.save_table_button.setEnabled(False)
        button_layout.addWidget(self.save_table_button)
        button_layout.addStretch()
        self.table_widget = QTableWidget()
        self.table_widget.setAlternatingRowColors(True)
        self.table_widget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        layout.addLayout(button_layout)
        layout.addWidget(self.table_widget)

    def init_manual_tab(self):
        layout = QVBoxLayout(self.tab4)
        manual_text = QTextEdit()
        manual_text.setReadOnly(True)
        manual_text.setFont(QFont("Segoe UI", 10))
        manual_content = """
        <h1>Steel Tempering Curve Analysis - User Manual</h1>
        <h2>Program Description</h2>
        <p>This program is designed for visualization and analysis of steel tempering process data.</p>
        <h2>License Information</h2>
        <p><b>GNU GENERAL PUBLIC LICENSE Version 3</b></p>
        <h2>Data File Requirements</h2>
        <ul><li>Steel type, Tempering time (s), Tempering temperature (ºC), Final hardness (HRC) - post tempering</li></ul>
        <p>Copyright (c) 2025-2026 Marek Góral m_goral@interia.pl</p>
        """
        manual_text.setHtml(manual_content)
        layout.addWidget(manual_text)

    def initialize_plot(self):
        self.ax.clear()
        self.ax.set_xlabel('Tempering Temperature [°C]')
        self.ax.set_ylabel('Final Hardness [HRC]')
        self.ax.set_title('Steel Tempering Curves')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        self.canvas.draw()

    def try_auto_load_data(self):
        default_file = "Tempering data for carbon and low alloy steels - Raiipa.csv"
        if os.path.exists(default_file):
            try:
                self.data = pd.read_csv(default_file)
                required_columns = ['Steel type', 'Tempering temperature (ºC)', 'Final hardness (HRC) - post tempering']
                if not all(col in self.data.columns for col in required_columns):
                    QMessageBox.critical(self, "Error", "CSV file does not contain required columns!")
                    return
                steel_types = sorted(self.data['Steel type'].unique())
                self.steel_combo.clear()
                self.steel_combo.addItems(steel_types)
                if steel_types:
                    self.steel_combo.setCurrentText(steel_types[0])
                    self.save_plot_button.setEnabled(True)
                    self.save_table_button.setEnabled(True)
                QMessageBox.information(self, "Success", f"Data loaded successfully from:\n{default_file}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load data: {str(e)}")
        else:
            QMessageBox.information(self, "Information", "Default data file not found. Use 'Load CSV Data' button.")

    def load_data(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select CSV File", "", "CSV Files (*.csv)")
        if file_path:
            try:
                self.data = pd.read_csv(file_path)
                required_columns = ['Steel type', 'Tempering temperature (ºC)', 'Final hardness (HRC) - post tempering']
                if not all(col in self.data.columns for col in required_columns):
                    QMessageBox.critical(self, "Error", "CSV file does not contain required columns!")
                    return
                steel_types = sorted(self.data['Steel type'].unique())
                self.steel_combo.clear()
                self.steel_combo.addItems(steel_types)
                if steel_types:
                    self.steel_combo.setCurrentText(steel_types[0])
                    self.save_plot_button.setEnabled(True)
                    self.save_table_button.setEnabled(True)
                QMessageBox.information(self, "Success", f"Data loaded successfully!\nNumber of rows: {len(self.data)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load data: {str(e)}")

    def update_all(self, steel_type):
        if self.data is None or steel_type == "":
            return
        self.current_steel_type = steel_type
        self.current_steel_data = self.data[self.data['Steel type'] == steel_type]
        if self.current_steel_data.empty:
            return
        self.update_plot()
        self.update_info()
        self.update_table()

    def update_plot(self):
        if self.current_steel_data is None:
            return
        steel_data = self.current_steel_data.sort_values('Tempering temperature (ºC)')
        tempering_times = steel_data['Tempering time (s)'].unique()
        self.ax.clear()
        self.plots = []
        self.plot_data = []
        for time in tempering_times:
            time_data = steel_data[steel_data['Tempering time (s)'] == time]
            plot = self.ax.plot(time_data['Tempering temperature (ºC)'], 
                                time_data['Final hardness (HRC) - post tempering'],
                                marker='o', label=f'{time} s')
            self.plots.append(plot[0])
            for i, (temp, hardness) in enumerate(zip(time_data['Tempering temperature (ºC)'], 
                                                     time_data['Final hardness (HRC) - post tempering'])):
                self.plot_data.append({'x': temp, 'y': hardness, 'time': time, 'plot': plot[0]})
        self.ax.set_xlabel('Tempering Temperature [°C]')
        self.ax.set_ylabel('Final Hardness [HRC]')
        self.ax.set_title(f'Tempering Curves for {self.current_steel_type}')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        self.ax.legend(title='Tempering Time')
        self.ax.autoscale()
        self.canvas.draw()

    def update_info(self):
        if self.current_steel_data is None:
            return
        first_row = self.current_steel_data.iloc[0]
        info_text = f"<h2>Steel Information: {self.current_steel_type}</h2><h3>Chemical Composition:</h3><ul>"
        elements = ['C (%wt)', 'Mn (%wt)', 'P (%wt)', 'S (%wt)', 'Si (%wt)', 
                    'Ni (%wt)', 'Cr (%wt)', 'Mo (%wt)', 'V (%wt)', 'Al (%wt)', 'Cu (%wt)']
        for element in elements:
            value = first_row.get(element, 'no data')
            if value != 'no data' and not pd.isna(value):
                element_name = element.split(' ')[0]
                info_text += f"<li>{element_name}: {value}%</li>"
        info_text += "</ul><h3>Additional Information:</h3>"
        info_text += f"<p>Number of measurements: {len(self.current_steel_data)}</p>"
        min_temp = self.current_steel_data['Tempering temperature (ºC)'].min()
        max_temp = self.current_steel_data['Tempering temperature (ºC)'].max()
        min_hard = self.current_steel_data['Final hardness (HRC) - post tempering'].min()
        max_hard = self.current_steel_data['Final hardness (HRC) - post tempering'].max()
        info_text += f"<p>Temperature range: {min_temp} - {max_temp}°C</p>"
        info_text += f"<p>Hardness range: {min_hard} - {max_hard} HRC</p>"
        self.info_text.setHtml(info_text)

    def update_table(self):
        if self.current_steel_data is None:
            return
        sorted_data = self.current_steel_data.sort_values(['Tempering time (s)', 'Tempering temperature (ºC)'])
        self.table_widget.setRowCount(len(sorted_data))
        self.table_widget.setColumnCount(4)
        self.table_widget.setHorizontalHeaderLabels(['Time [s]', 'Temperature [°C]', 'Hardness [HRC]', 'Source'])
        for row_idx, (_, row_data) in enumerate(sorted_data.iterrows()):
            self.table_widget.setItem(row_idx, 0, QTableWidgetItem(str(row_data['Tempering time (s)'])))
            self.table_widget.setItem(row_idx, 1, QTableWidgetItem(str(row_data['Tempering temperature (ºC)'])))
            self.table_widget.setItem(row_idx, 2, QTableWidgetItem(str(row_data['Final hardness (HRC) - post tempering'])))
            self.table_widget.setItem(row_idx, 3, QTableWidgetItem(str(row_data.get('Source', 'no data'))))

    def on_hover(self, event):
        if event.inaxes != self.ax or not self.plot_data:
            if self.annotation:
                self.annotation.remove()
                self.annotation = None
                self.canvas.draw()
            return
        min_dist = float('inf')
        closest = None
        for point in self.plot_data:
            dist = ((event.xdata - point['x'])**2 + (event.ydata - point['y'])**2)**0.5
            if dist < min_dist:
                min_dist = dist
                closest = point
        if min_dist < 15:
            if self.annotation:
                self.annotation.remove()
            self.annotation = self.ax.annotate(
                f'Time: {closest["time"]} s\nTemp: {closest["x"]:.1f}°C\nHardness: {closest["y"]:.1f} HRC',
                xy=(closest["x"], closest["y"]),
                xytext=(20, 20), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', alpha=0.9),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
            )
            self.canvas.draw()
        elif self.annotation:
            self.annotation.remove()
            self.annotation = None
            self.canvas.draw()

    def save_plot(self):
        if self.current_steel_data is None:
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Plot", f"plot_{self.current_steel_type}.png",
                                                   "Images (*.png *.jpg *.pdf *.svg)")
        if file_path:
            try:
                self.figure.savefig(file_path, dpi=300, bbox_inches='tight')
                QMessageBox.information(self, "Success", f"Plot saved to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed: {str(e)}")

    def save_table(self):
        if self.current_steel_data is None:
            return
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Table", f"table_{self.current_steel_type}.csv",
                                                   "CSV Files (*.csv)")
        if file_path:
            try:
                sorted_data = self.current_steel_data.sort_values(['Tempering time (s)', 'Tempering temperature (ºC)'])
                sorted_data.to_csv(file_path, index=False, encoding='utf-8')
                QMessageBox.information(self, "Success", f"Table saved to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed: {str(e)}")

# =============================================================================
# 2. JOMINA ANALYSER – Hardenability calculations (zwiększona czcionka, ikony)
# =============================================================================
STANDARD_JOMINY_DISTANCES = [1.5, 3, 5, 7, 9, 11, 13, 15, 20, 25, 30, 35, 40, 45, 50]

def hrc_to_hb(hrc):
    if hrc < 20:
        return hrc * 10
    elif hrc > 60:
        return 650
    else:
        return 0.363 * hrc**2 + 13.93 * hrc + 100.8

def hb_to_hrc(hb):
    if hb < 200:
        return hb / 10
    elif hb > 600:
        return 65
    else:
        return 0.094 * hb - 16.9

def half_martensite_hardness(carbon_content):
    if carbon_content < 0.1:
        return 26
    elif carbon_content < 0.11:
        return 27
    elif carbon_content < 0.12:
        return 27
    elif carbon_content < 0.13:
        return 28
    elif carbon_content < 0.14:
        return 28
    elif carbon_content < 0.15:
        return 29
    elif carbon_content < 0.16:
        return 30
    elif carbon_content < 0.17:
        return 30
    elif carbon_content < 0.18:
        return 31
    elif carbon_content < 0.19:
        return 31
    elif carbon_content < 0.20:
        return 32
    elif carbon_content < 0.21:
        return 32
    elif carbon_content < 0.22:
        return 33
    elif carbon_content < 0.23:
        return 34
    elif carbon_content < 0.24:
        return 34
    elif carbon_content < 0.25:
        return 35
    elif carbon_content < 0.26:
        return 35
    elif carbon_content < 0.27:
        return 36
    elif carbon_content < 0.28:
        return 36
    elif carbon_content < 0.29:
        return 37
    elif carbon_content < 0.30:
        return 37
    elif carbon_content < 0.31:
        return 38
    elif carbon_content < 0.32:
        return 38
    elif carbon_content < 0.33:
        return 39
    elif carbon_content < 0.34:
        return 40
    elif carbon_content < 0.35:
        return 40
    elif carbon_content < 0.36:
        return 41
    elif carbon_content < 0.37:
        return 41
    elif carbon_content < 0.38:
        return 42
    elif carbon_content < 0.39:
        return 42
    elif carbon_content < 0.40:
        return 43
    elif carbon_content < 0.41:
        return 43
    elif carbon_content < 0.42:
        return 43
    elif carbon_content < 0.43:
        return 44
    elif carbon_content < 0.44:
        return 44
    elif carbon_content < 0.45:
        return 45
    elif carbon_content < 0.46:
        return 45
    elif carbon_content < 0.47:
        return 45
    elif carbon_content < 0.48:
        return 46
    elif carbon_content < 0.49:
        return 46
    elif carbon_content < 0.50:
        return 47
    elif carbon_content < 0.51:
        return 47
    elif carbon_content < 0.52:
        return 48
    elif carbon_content < 0.53:
        return 48
    elif carbon_content < 0.54:
        return 48
    elif carbon_content < 0.55:
        return 49
    elif carbon_content < 0.56:
        return 49
    elif carbon_content < 0.57:
        return 50
    elif carbon_content < 0.58:
        return 50
    elif carbon_content < 0.59:
        return 51
    elif carbon_content < 0.60:
        return 51
    elif carbon_content < 0.61:
        return 51
    elif carbon_content < 0.62:
        return 51
    elif carbon_content < 0.63:
        return 52
    elif carbon_content < 0.64:
        return 52
    elif carbon_content < 0.65:
        return 52
    elif carbon_content < 0.66:
        return 52
    elif carbon_content < 0.67:
        return 53
    elif carbon_content < 0.68:
        return 53
    elif carbon_content < 0.69:
        return 53
    else:
        return 53

def astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu=0):
    if C < 0.1:
        D_I_base = 0.05
    elif C < 0.2:
        D_I_base = 0.06
    elif C < 0.3:
        D_I_base = 0.08
    elif C < 0.4:
        D_I_base = 0.13
    elif C < 0.5:
        D_I_base = 0.18
    elif C < 0.6:
        D_I_base = 0.22
    elif C < 0.7:
        D_I_base = 0.26
    elif C < 0.8:
        D_I_base = 0.30
    elif C < 0.9:
        D_I_base = 0.34
    else:
        D_I_base = 0.38
    f_Mn = 1.0 + 4.10 * Mn
    f_Si = 1.0 + 0.70 * Si
    f_Ni = 1.0 + 0.50 * Ni
    f_Cr = 1.0 + 2.33 * Cr
    f_Mo = 1.0 + 3.14 * Mo
    f_Cu = 1.0 + 0.52 * Cu
    DI = D_I_base * f_Mn * f_Si * f_Ni * f_Cr * f_Mo * f_Cu
    return DI * 25.4

def just_DI(C, Mn, Cr, Mo, Ni, V):
    K = 25.4
    f_C = 0.123 + 0.289*C + 0.330*C**2
    f_Mn = 1 + 3.333*Mn
    f_Cr = 1 + 1.167*Cr
    f_Mo = 1 + 2.5*Mo
    f_Ni = 1 + 0.362*Ni
    f_V = 1 + 0.667*V
    return K * f_C * f_Mn * f_Cr * f_Mo * f_Ni * f_V

def cremona_DI(C, Mn, Si, Cr, Ni, Mo, Cu, V, W):
    K = 25.4
    f_C = 0.19 + 0.27*C
    f_Mn = 1 + 3*Mn
    f_Si = 1 + 0.7*Si
    f_Cr = 1 + Cr
    f_Ni = 1 + 0.27*Ni
    f_Mo = 1 + 3*Mo
    f_Cu = 1 + 0.4*Cu
    f_V = 1 + 2*V
    f_W = 1 + 2*W
    return K * f_C * f_Mn * f_Si * f_Cr * f_Ni * f_Mo * f_Cu * f_V * f_W

def just_jominy_curve(C, Mn, Si, Cr, Ni, Mo, V, grain_size, distances):
    m = 1 / 1.5875
    hardness = []
    for e in distances:
        term1 = 88 * C
        term2 = -0.0135 * (m * e)**2 / max(C, 0.001)
        term3 = 19 * Cr
        term4 = 6.3 * Ni
        term5 = 16 * Mn
        term6 = 35 * Mo
        term7 = 5 * Si
        term8 = -0.82 * grain_size
        term9 = -20 * np.sqrt(m * e)
        term10 = 2.11 * (m * e)
        constant = -2
        hv = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + constant
        hardness.append(max(0, min(65, hv)))
    return np.array(hardness)

def astm_jominy_curve(C, Mn, Si, Ni, Cr, Mo, Cu, grain_size, distances):
    DI = astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu)
    if C < 0.10:
        initial_hardness = 38
    elif C < 0.11:
        initial_hardness = 39
    elif C < 0.12:
        initial_hardness = 40
    elif C < 0.13:
        initial_hardness = 40
    elif C < 0.14:
        initial_hardness = 41
    elif C < 0.15:
        initial_hardness = 41
    elif C < 0.16:
        initial_hardness = 42
    elif C < 0.17:
        initial_hardness = 43
    elif C < 0.18:
        initial_hardness = 43
    elif C < 0.19:
        initial_hardness = 44
    elif C < 0.20:
        initial_hardness = 44
    elif C < 0.21:
        initial_hardness = 45
    elif C < 0.22:
        initial_hardness = 46
    elif C < 0.23:
        initial_hardness = 46
    elif C < 0.24:
        initial_hardness = 46
    elif C < 0.25:
        initial_hardness = 47
    elif C < 0.26:
        initial_hardness = 48
    elif C < 0.27:
        initial_hardness = 49
    elif C < 0.28:
        initial_hardness = 49
    elif C < 0.29:
        initial_hardness = 50
    elif C < 0.30:
        initial_hardness = 50
    elif C < 0.31:
        initial_hardness = 51
    elif C < 0.32:
        initial_hardness = 51
    elif C < 0.33:
        initial_hardness = 52
    elif C < 0.34:
        initial_hardness = 53
    elif C < 0.35:
        initial_hardness = 53
    elif C < 0.36:
        initial_hardness = 54
    elif C < 0.37:
        initial_hardness = 55
    elif C < 0.38:
        initial_hardness = 55
    elif C < 0.39:
        initial_hardness = 56
    elif C < 0.40:
        initial_hardness = 56
    elif C < 0.41:
        initial_hardness = 57
    elif C < 0.42:
        initial_hardness = 57
    elif C < 0.43:
        initial_hardness = 58
    elif C < 0.44:
        initial_hardness = 58
    elif C < 0.45:
        initial_hardness = 59
    elif C < 0.46:
        initial_hardness = 59
    elif C < 0.47:
        initial_hardness = 59
    elif C < 0.48:
        initial_hardness = 59
    elif C < 0.49:
        initial_hardness = 60
    elif C < 0.50:
        initial_hardness = 61
    elif C < 0.51:
        initial_hardness = 61
    elif C < 0.52:
        initial_hardness = 62
    elif C < 0.53:
        initial_hardness = 62
    elif C < 0.54:
        initial_hardness = 63
    elif C < 0.55:
        initial_hardness = 63
    elif C < 0.56:
        initial_hardness = 63
    elif C < 0.57:
        initial_hardness = 64
    elif C < 0.58:
        initial_hardness = 64
    elif C < 0.59:
        initial_hardness = 64
    elif C < 0.60:
        initial_hardness = 64
    elif C < 0.61:
        initial_hardness = 64
    elif C < 0.62:
        initial_hardness = 65
    elif C < 0.63:
        initial_hardness = 65
    elif C < 0.64:
        initial_hardness = 65
    elif C < 0.65:
        initial_hardness = 65
    elif C < 0.66:
        initial_hardness = 65
    elif C < 0.67:
        initial_hardness = 65
    elif C < 0.68:
        initial_hardness = 65
    elif C < 0.69:
        initial_hardness = 65
    else:
        initial_hardness = 65
    half_martensitic = half_martensite_hardness(C)
    hardness = []
    for dist in distances:
        decay = np.exp(-dist / (DI / 10))
        hv = half_martensitic + (initial_hardness - half_martensitic) * decay
        hardness.append(max(20, min(65, hv)))
    return np.array(hardness)

def calculate_ideal_DI(actual_DI, severity_factor=1.0):
    return actual_DI * (1.0 / 0.8) ** 0.5

def jominy_distance_to_DI(j_distance):
    dist_data = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,
                 16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,
                 31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,
                 46.0,47.0,48.0,49.0,50.0]
    di_data = [8.4,15.7,22.9,29.7,36.3,42.9,48.2,54.2,59.5,64.2,68.6,72.1,76.4,80.1,84.0,
               87.6,90.1,94.2,97.1,100.5,103.7,106.5,109.7,112.2,114.9,117.4,119.9,122.4,124.7,127.1,
               129.0,131.4,133.5,135.2,137.1,139.1,140.9,142.8,144.7,146.4,148.3,150.1,151.7,153.4,154.1,
               156.5,157.8,159.2,160.5,161.8]
    return np.interp(j_distance, dist_data, di_data)

steel_database_jomina = {
    "1.6580; 30CrNiMo8; 30H2N2M": {
        "C": 0.30, "Mn": 0.50, "Si": 0.30, "Cr": 2.00, "Ni": 2.00, "Mo": 0.30, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 120, "half_martensitic": 37,
        "description": "Heat treatable steel, used for shafts, axles, gears, high-strength bolts",
        "heat_treatment": "Quenching: 850-880°C in oil, Tempering: 550-650°C",
        "hardness": "280-320 HB",
        "hardness_after_quenching": "50-55 HRC"
    },
    "1.6582; 34CrNiMo6; 34HNM": {
        "C": 0.34, "Mn": 0.60, "Si": 0.30, "Cr": 1.50, "Ni": 1.50, "Mo": 0.25, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 100, "half_martensitic": 40,
        "description": "Heat treatable steel, used for high-strength machine components, crankshafts",
        "heat_treatment": "Quenching: 840-860°C in oil, Tempering: 580-650°C",
        "hardness": "300-340 HB",
        "hardness_after_quenching": "52-56 HRC"
    },
    "1.7035; 41Cr4; 40H": {
        "C": 0.41, "Mn": 0.70, "Si": 0.30, "Cr": 1.00, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 80, "half_martensitic": 43,
        "description": "Case hardening and heat treatable steel",
        "heat_treatment": "Quenching: 830-850°C in oil, Tempering: 550-650°C",
        "hardness": "200-250 HB",
        "hardness_after_quenching": "45-50 HRC"
    },
    "1.7225; 42CrMo4; 40HM": {
        "C": 0.42, "Mn": 0.75, "Si": 0.30, "Cr": 1.00, "Ni": 0.10, "Mo": 0.20, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 90, "half_martensitic": 43,
        "description": "Heat treatable steel, used for shafts, axles, fasteners, gears",
        "heat_treatment": "Quenching: 840-860°C in oil, Tempering: 550-650°C",
        "hardness": "280-320 HB",
        "hardness_after_quenching": "48-52 HRC"
    },
    "1.0503; C45; 45": {
        "C": 0.45, "Mn": 0.65, "Si": 0.25, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 30, "half_martensitic": 45,
        "description": "General purpose structural steel",
        "heat_treatment": "Quenching: 820-850°C in water, Tempering: 550-600°C",
        "hardness": "170-220 HB",
        "hardness_after_quenching": "25-30 HRC"
    },
    "1.0535; C55; 55": {
        "C": 0.55, "Mn": 0.65, "Si": 0.25, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 35, "half_martensitic": 49,
        "description": "Higher strength structural steel",
        "heat_treatment": "Quenching: 820-850°C in water, Tempering: 450-550°C",
        "hardness": "200-250 HB",
        "hardness_after_quenching": "30-35 HRC"
    },
    "1.0545; S355; 18G2A": {
        "C": 0.20, "Mn": 1.40, "Si": 0.40, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 25, "half_martensitic": 32,
        "description": "General purpose structural steel",
        "heat_treatment": "Normalizing: 890-950°C",
        "hardness": "140-180 HB",
        "hardness_after_quenching": "Not applicable"
    },
    "1.0601; C60; 60": {
        "C": 0.60, "Mn": 0.75, "Si": 0.25, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 40, "half_martensitic": 51,
        "description": "Higher strength structural steel, used for springs",
        "heat_treatment": "Quenching: 810-840°C in water, Tempering: 400-500°C",
        "hardness": "220-270 HB",
        "hardness_after_quenching": "35-40 HRC"
    },
    "1.6587; 18CrNiMo7-6; 17HNM": {
        "C": 0.18, "Mn": 0.60, "Si": 0.30, "Cr": 1.70, "Ni": 1.50, "Mo": 0.30, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 110, "half_martensitic": 31,
        "description": "Case hardening steel, used for gears, shafts",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (surface), 280-320 HB (core)",
        "hardness_after_quenching": "58-62 HRC"
    },
    "1.7131; 16MnCr5; 16HG": {
        "C": 0.16, "Mn": 1.10, "Si": 0.30, "Cr": 1.00, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 70, "half_martensitic": 30,
        "description": "Case hardening steel, used for gears",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (surface), 250-300 HB (core)",
        "hardness_after_quenching": "58-62 HRC"
    },
    "1.7147; 20MnCr5; 20HG": {
        "C": 0.20, "Mn": 1.10, "Si": 0.30, "Cr": 1.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 75, "half_martensitic": 32,
        "description": "Case hardening steel",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (surface), 280-320 HB (core)",
        "hardness_after_quenching": "58-62 HRC"
    },
    "1.8509; 41CrAlMo7-10; 38HMJ": {
        "C": 0.41, "Mn": 0.60, "Si": 0.30, "Cr": 1.50, "Ni": 0.10, "Mo": 0.25, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 85, "half_martensitic": 43,
        "description": "Nitriding steel",
        "heat_treatment": "Quenching: 930-950°C in oil, Tempering: 650-700°C",
        "hardness": "280-320 HB",
        "hardness_after_quenching": "45-50 HRC"
    },
}

class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About Jomina Analyzer")
        self.setModal(True)
        self.resize(500,400)
        layout = QVBoxLayout()
        title = QLabel("Jomina Analyzer")
        title.setStyleSheet("font-size: 18pt; font-weight: bold;")
        layout.addWidget(title)
        version = QLabel("Version 1.0")
        layout.addWidget(version)
        desc = QTextEdit()
        desc.setReadOnly(True)
        desc.setPlainText("Jominy Analyzer - Hardenability calculations using ASTM/Grossman, Just, de Cremona methods.")
        layout.addWidget(desc)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
        self.setLayout(layout)

class JominyPlot(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        self.axes = self.fig.add_subplot(111)
        self.axes.grid(True)
        self.axes.set_xlabel('Distance from quenched end [mm]')
        self.axes.set_ylabel('Hardness [HRC]')
        self.axes.set_title('Jominy Curve')
        self.axes.set_xlim(0,50)
        self.curves = []

    def plot_jominy(self, dist_just, hard_just, dist_astm, hard_astm, experimental_data, steel_name, half_mart, show_just):
        self.curves = []
        if show_just:
            self.curves.append({'distance': dist_just, 'hardness': hard_just, 'curve_name': "Just Method", 'color':'r', 'linestyle':'-'})
        self.curves.append({'distance': dist_astm, 'hardness': hard_astm, 'curve_name': "ASTM Method", 'color':'b', 'linestyle':'--'})
        if experimental_data:
            colors = ['g','m','c','y']
            markers = ['o','s','^','D']
            for i, (d,h) in enumerate(experimental_data):
                if len(d)>0:
                    self.curves.append({'distance': d, 'hardness': h, 'curve_name': f"Exp {i+1}", 'color':colors[i%4], 'linestyle':'-', 'marker':markers[i%4]})
        self.update_plot(half_mart)

    def update_plot(self, half_mart=None):
        self.axes.clear()
        for curve in self.curves:
            if 'marker' in curve:
                self.axes.plot(curve['distance'], curve['hardness'], color=curve['color'], linestyle=curve['linestyle'],
                               marker=curve['marker'], markersize=4, linewidth=2, label=curve['curve_name'])
            else:
                self.axes.plot(curve['distance'], curve['hardness'], color=curve['color'], linestyle=curve['linestyle'],
                               linewidth=2, label=curve['curve_name'])
        if half_mart:
            self.axes.axhline(y=half_mart, color='g', linestyle='--', label=f'Half-martensitic: {half_mart} HRC')
        self.axes.grid(True)
        self.axes.set_xlabel('Distance [mm]')
        self.axes.set_ylabel('Hardness [HRC]')
        self.axes.set_xlim(0,50)
        self.axes.legend()
        self.draw()

class DIComparisonPlot(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        self.axes = self.fig.add_subplot(111)

    def plot_comparison(self, methods, di_vals, steel_name, jominy_di, ideal_di):
        self.axes.clear()
        ypos = np.arange(len(methods))
        self.axes.barh(ypos, di_vals, align='center', alpha=0.7)
        self.axes.set_yticks(ypos)
        self.axes.set_yticklabels(methods)
        self.axes.set_xlabel('Critical Diameter DI [mm]')
        self.axes.set_title(f'DI Comparison for {steel_name}')
        for i,v in enumerate(di_vals):
            self.axes.text(v+1, i, f'{v:.2f}', va='center')
        if jominy_di:
            self.axes.axvline(x=jominy_di, color='r', linestyle='--', label=f'Actual DI: {jominy_di:.2f} mm')
        if ideal_di:
            self.axes.axvline(x=ideal_di, color='g', linestyle='--', label=f'Ideal DI: {ideal_di:.2f} mm')
        if jominy_di or ideal_di:
            self.axes.legend()
        self.draw()

class SteelCalculatorApp(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Steel Hardenability Calculator - Jominy Analyzer")
        self.setGeometry(100,100,1600,900)
        # Zwiększenie czcionki
        font = QFont("Segoe UI", 10)
        self.setFont(font)
        
        self.jominy_curves = []
        self.jominy_di = None
        self.ideal_di = None
        self.show_just_method = True
        self.di_results = {}

        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tab4 = QWidget()
        self.tab5 = QWidget()
        self.tabs.addTab(self.tab1, "DI Calculations")
        self.tabs.addTab(self.tab2, "Jominy Curve")
        self.tabs.addTab(self.tab3, "DI Methods Comparison")
        self.tabs.addTab(self.tab4, "Methods Description")
        self.tabs.addTab(self.tab5, "About")
        self.setup_tab1()
        self.setup_tab2()
        self.setup_tab3()
        self.setup_tab4()
        self.setup_tab5()
        main_layout.addWidget(self.tabs)
        self.on_steel_changed(self.steel_combo.currentText())

    def setup_tab1(self):
        layout = QHBoxLayout(self.tab1)
        left = QVBoxLayout()
        steel_group = QGroupBox("Steel Selection")
        steel_layout = QVBoxLayout(steel_group)
        self.steel_combo = QComboBox()
        self.steel_combo.addItems(steel_database_jomina.keys())
        self.steel_combo.currentTextChanged.connect(self.on_steel_changed)
        steel_layout.addWidget(QLabel("Select steel:"))
        steel_layout.addWidget(self.steel_combo)
        left.addWidget(steel_group)

        comp_group = QGroupBox("Chemical Composition [%]")
        comp_layout = QFormLayout(comp_group)
        self.c_input = QLineEdit("0.30")
        self.mn_input = QLineEdit("0.50")
        self.si_input = QLineEdit("0.30")
        self.cr_input = QLineEdit("2.00")
        self.ni_input = QLineEdit("2.00")
        self.mo_input = QLineEdit("0.30")
        self.cu_input = QLineEdit("0.10")
        self.v_input = QLineEdit("0.00")
        self.w_input = QLineEdit("0.00")
        comp_layout.addRow("C:", self.c_input)
        comp_layout.addRow("Mn:", self.mn_input)
        comp_layout.addRow("Si:", self.si_input)
        comp_layout.addRow("Cr:", self.cr_input)
        comp_layout.addRow("Ni:", self.ni_input)
        comp_layout.addRow("Mo:", self.mo_input)
        comp_layout.addRow("Cu:", self.cu_input)
        comp_layout.addRow("V:", self.v_input)
        comp_layout.addRow("W:", self.w_input)
        left.addWidget(comp_group)

        params_group = QGroupBox("Parameters")
        params_layout = QFormLayout(params_group)
        self.grain_size_input = QLineEdit("7")
        self.hardness_after_quenching_input = QLineEdit("50")
        params_layout.addRow("ASTM Grain Size:", self.grain_size_input)
        params_layout.addRow("Hardness after quenching [HRC]:", self.hardness_after_quenching_input)
        left.addWidget(params_group)

        btn_layout = QHBoxLayout()
        self.calc_button = QPushButton(" Calculate DI")
        self.calc_button.setIcon(self.style().standardIcon(QStyle.SP_ComputerIcon))
        self.calc_button.clicked.connect(self.calculate_DI)
        self.save_comp_button = QPushButton(" Save Composition")
        self.save_comp_button.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_comp_button.clicked.connect(self.save_composition)
        self.load_comp_button = QPushButton(" Load Composition")
        self.load_comp_button.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.load_comp_button.clicked.connect(self.load_composition)
        btn_layout.addWidget(self.calc_button)
        btn_layout.addWidget(self.save_comp_button)
        btn_layout.addWidget(self.load_comp_button)
        left.addLayout(btn_layout)
        left.addStretch()

        middle = QVBoxLayout()
        info_group = QGroupBox("Steel Information")
        info_layout = QVBoxLayout(info_group)
        self.steel_description = QTextEdit()
        self.steel_description.setReadOnly(True)
        info_layout.addWidget(self.steel_description)
        middle.addWidget(info_group)
        heat_group = QGroupBox("Heat Treatment")
        heat_layout = QVBoxLayout(heat_group)
        self.heat_treatment_info = QTextEdit()
        self.heat_treatment_info.setReadOnly(True)
        heat_layout.addWidget(self.heat_treatment_info)
        middle.addWidget(heat_group)
        middle.addStretch()

        right = QVBoxLayout()
        results_group = QGroupBox("DI Results [mm]")
        res_layout = QFormLayout(results_group)
        self.astm_grossman_result = QLabel("—")
        self.just_result = QLabel("—")
        self.cremona_result = QLabel("—")
        self.experimental_result = QLabel("—")
        self.ideal_experimental_result = QLabel("—")
        res_layout.addRow("ASTM/Grossman:", self.astm_grossman_result)
        res_layout.addRow("Just:", self.just_result)
        res_layout.addRow("de Cremona:", self.cremona_result)
        res_layout.addRow("Experimental (actual):", self.experimental_result)
        res_layout.addRow("Experimental (ideal):", self.ideal_experimental_result)
        right.addWidget(results_group)
        self.save_results_button = QPushButton(" Save Results")
        self.save_results_button.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_results_button.clicked.connect(self.save_results)
        right.addWidget(self.save_results_button)
        right.addStretch()

        left_w = QWidget(); left_w.setLayout(left); left_w.setMaximumWidth(400)
        mid_w = QWidget(); mid_w.setLayout(middle); mid_w.setMaximumWidth(350)
        right_w = QWidget(); right_w.setLayout(right); right_w.setMaximumWidth(300)
        layout.addWidget(left_w); layout.addWidget(mid_w); layout.addWidget(right_w)

    def setup_tab2(self):
        layout = QHBoxLayout(self.tab2)
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        jominy_group = QGroupBox("Jominy Curve Data")
        jominy_layout = QVBoxLayout(jominy_group)
        curve_select = QHBoxLayout()
        self.curve_checkboxes = []
        for i in range(4):
            cb = QCheckBox(f"Curve {i+1}")
            cb.setChecked(i==0)
            cb.toggled.connect(self.update_jominy_plot)
            self.curve_checkboxes.append(cb)
            curve_select.addWidget(cb)
        self.show_just_checkbox = QCheckBox("Show Just Method")
        self.show_just_checkbox.setChecked(True)
        self.show_just_checkbox.toggled.connect(self.toggle_just_method)
        curve_select.addWidget(self.show_just_checkbox)
        jominy_layout.addLayout(curve_select)

        self.jominy_table = QTableWidget()
        self.jominy_table.setColumnCount(5)
        self.jominy_table.setHorizontalHeaderLabels(["Distance [mm]", "Hardness 1", "Hardness 2", "Hardness 3", "Hardness 4"])
        self.jominy_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.jominy_table.setMinimumHeight(400)

        btn_row = QHBoxLayout()
        self.add_row_button = QPushButton(" Add Row")
        self.add_row_button.setIcon(self.style().standardIcon(QStyle.SP_FileDialogNewFolder))
        self.add_row_button.clicked.connect(self.add_jominy_row)
        self.remove_row_button = QPushButton(" Remove Row")
        self.remove_row_button.setIcon(self.style().standardIcon(QStyle.SP_TrashIcon))
        self.remove_row_button.clicked.connect(self.remove_jominy_row)
        self.calc_di_button = QPushButton(" Calculate DI from Curve")
        self.calc_di_button.setIcon(self.style().standardIcon(QStyle.SP_ComputerIcon))
        self.calc_di_button.clicked.connect(self.calculate_DI_from_jominy)
        btn_row.addWidget(self.add_row_button)
        btn_row.addWidget(self.remove_row_button)
        btn_row.addWidget(self.calc_di_button)
        jominy_layout.addWidget(self.jominy_table)
        jominy_layout.addLayout(btn_row)

        file_group = QGroupBox("File Operations")
        file_layout = QHBoxLayout(file_group)
        self.load_jominy_button = QPushButton(" Load Jominy Data")
        self.load_jominy_button.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.load_jominy_button.clicked.connect(self.load_jominy_data)
        self.save_jominy_button = QPushButton(" Save Jominy Data")
        self.save_jominy_button.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_jominy_button.clicked.connect(self.save_jominy_data)
        self.save_plot_button = QPushButton(" Save Plot")
        self.save_plot_button.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_plot_button.clicked.connect(self.save_jominy_plot)
        file_layout.addWidget(self.load_jominy_button)
        file_layout.addWidget(self.save_jominy_button)
        file_layout.addWidget(self.save_plot_button)
        left_layout.addWidget(jominy_group)
        left_layout.addWidget(file_group)
        left_layout.addStretch()

        self.jominy_plot = JominyPlot(self, width=8, height=8, dpi=100)
        layout.addWidget(left_panel)
        layout.addWidget(self.jominy_plot)
        self.fill_standard_distances()
        self.jominy_table.cellChanged.connect(self.update_jominy_plot)

    def fill_standard_distances(self):
        self.jominy_table.setRowCount(len(STANDARD_JOMINY_DISTANCES))
        for i, d in enumerate(STANDARD_JOMINY_DISTANCES):
            self.jominy_table.setItem(i,0, QTableWidgetItem(str(d)))
            for j in range(1,5):
                self.jominy_table.setItem(i,j, QTableWidgetItem(""))

    def setup_tab3(self):
        layout = QVBoxLayout(self.tab3)
        self.di_comparison_plot = DIComparisonPlot(self, width=12, height=8, dpi=100)
        layout.addWidget(self.di_comparison_plot)
        self.update_comparison_button = QPushButton(" Update Comparison")
        self.update_comparison_button.setIcon(self.style().standardIcon(QStyle.SP_BrowserReload))
        self.update_comparison_button.clicked.connect(self.update_di_comparison)
        layout.addWidget(self.update_comparison_button)

    def setup_tab4(self):
        layout = QVBoxLayout(self.tab4)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText("""Methods description:
ASTM/Grossman: DI = base * f_Mn * f_Si * f_Ni * f_Cr * f_Mo * f_Cu
Just (1969): DI = K * f_C * f_Mn * f_Cr * f_Mo * f_Ni * f_V
de Cremona (1970): DI = K * f_C * f_Mn * f_Si * f_Cr * f_Ni * f_Mo * f_Cu * f_V * f_W
Jominy curve calculations: Just method and ASTM A255 method.
""")
        layout.addWidget(text)

    def setup_tab5(self):
        layout = QVBoxLayout(self.tab5)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setHtml("<h1>Jomina Analyzer 1.2</h1><p>GNU GPL v3</p><p>Contact: m_goral@interia.pl</p>")
        layout.addWidget(text)

    def toggle_just_method(self, checked):
        self.show_just_method = checked
        self.update_jominy_plot()

    def on_steel_changed(self, steel_name):
        if steel_name in steel_database_jomina:
            comp = steel_database_jomina[steel_name]
            self.c_input.setText(str(comp["C"]))
            self.mn_input.setText(str(comp["Mn"]))
            self.si_input.setText(str(comp["Si"]))
            self.cr_input.setText(str(comp["Cr"]))
            self.ni_input.setText(str(comp["Ni"]))
            self.mo_input.setText(str(comp["Mo"]))
            self.cu_input.setText(str(comp["Cu"]))
            self.v_input.setText(str(comp["V"]))
            self.w_input.setText(str(comp["W"]))
            self.grain_size_input.setText(str(comp.get("grain_size",7)))
            self.steel_description.setPlainText(comp.get("description",""))
            self.heat_treatment_info.setPlainText(comp.get("heat_treatment",""))
            self.update_jominy_plot()

    def get_composition(self):
        try:
            C = float(self.c_input.text())
            Mn = float(self.mn_input.text())
            Si = float(self.si_input.text())
            Cr = float(self.cr_input.text())
            Ni = float(self.ni_input.text())
            Mo = float(self.mo_input.text())
            Cu = float(self.cu_input.text())
            V = float(self.v_input.text())
            W = float(self.w_input.text())
            gs = float(self.grain_size_input.text())
            hq = float(self.hardness_after_quenching_input.text())
            return C, Mn, Si, Cr, Ni, Mo, Cu, V, W, gs, hq
        except:
            QMessageBox.warning(self, "Error", "Invalid input")
            return None

    def calculate_DI(self):
        comp = self.get_composition()
        if comp is None: return
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, gs, hq = comp
        di_astm = astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu)
        di_just = just_DI(C, Mn, Cr, Mo, Ni, V)
        di_crem = cremona_DI(C, Mn, Si, Cr, Ni, Mo, Cu, V, W)
        self.di_results = {"ASTM/Grossman": di_astm, "Just": di_just, "de Cremona": di_crem, "Steel": self.steel_combo.currentText()}
        self.astm_grossman_result.setText(f"{di_astm:.2f}")
        self.just_result.setText(f"{di_just:.2f}")
        self.cremona_result.setText(f"{di_crem:.2f}")
        self.update_di_comparison()
        self.update_jominy_plot()

    def update_di_comparison(self):
        comp = self.get_composition()
        if comp is None: return
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, gs, hq = comp
        di_astm = astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu)
        di_just = just_DI(C, Mn, Cr, Mo, Ni, V)
        di_crem = cremona_DI(C, Mn, Si, Cr, Ni, Mo, Cu, V, W)
        methods = ["ASTM/Grossman", "Just", "de Cremona"]
        values = [di_astm, di_just, di_crem]
        self.di_comparison_plot.plot_comparison(methods, values, self.steel_combo.currentText(), self.jominy_di, self.ideal_di)

    def add_jominy_row(self):
        row = self.jominy_table.rowCount()
        self.jominy_table.insertRow(row)
        for i in range(5):
            self.jominy_table.setItem(row, i, QTableWidgetItem(""))

    def remove_jominy_row(self):
        row = self.jominy_table.currentRow()
        if row>=0:
            self.jominy_table.removeRow(row)

    def get_jominy_data(self):
        distances = []
        hardness_curves = [[] for _ in range(4)]
        for row in range(self.jominy_table.rowCount()):
            dist_item = self.jominy_table.item(row,0)
            if dist_item and dist_item.text():
                try:
                    d = float(dist_item.text())
                    distances.append(d)
                    for i in range(4):
                        h_item = self.jominy_table.item(row, i+1)
                        if h_item and h_item.text():
                            hardness_curves[i].append(float(h_item.text()))
                        else:
                            hardness_curves[i].append(np.nan)
                except:
                    continue
        enabled = []
        for i in range(4):
            if self.curve_checkboxes[i].isChecked():
                valid = [ (distances[j], hardness_curves[i][j]) for j in range(len(distances)) if not np.isnan(hardness_curves[i][j]) ]
                if valid:
                    d_arr, h_arr = zip(*valid)
                    enabled.append((np.array(d_arr), np.array(h_arr)))
        return enabled

    def calculate_DI_from_jominy(self):
        curves = self.get_jominy_data()
        if not curves:
            QMessageBox.warning(self, "Error", "No experimental data")
            return
        comp = self.get_composition()
        if comp is None: return
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, gs, hq = comp
        dist, hard = curves[0]
        if len(dist)<2:
            QMessageBox.warning(self, "Insufficient data")
            return
        idx = np.argsort(dist)
        dist = dist[idx]
        hard = hard[idx]
        target = hq
        if hard[0] < target:
            di_val = 0
        else:
            below = np.where(hard < target)[0]
            if len(below)>0:
                first = below[0]
                if first==0:
                    di_val = 0
                else:
                    x1, x2 = dist[first-1], dist[first]
                    y1, y2 = hard[first-1], hard[first]
                    di_val = x1 + (x2-x1)*(target - y1)/(y2-y1)
            else:
                di_val = dist[-1]
        di_val = jominy_distance_to_DI(di_val)
        ideal = calculate_ideal_DI(di_val)
        self.jominy_di = di_val
        self.ideal_di = ideal
        self.experimental_result.setText(f"{di_val:.2f}")
        self.ideal_experimental_result.setText(f"{ideal:.2f}")
        self.update_jominy_plot()
        self.update_di_comparison()

    def update_jominy_plot(self):
        comp = self.get_composition()
        if comp is None: return
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, gs, hq = comp
        dist_just = np.linspace(0,50,100)
        hard_just = just_jominy_curve(C, Mn, Si, Cr, Ni, Mo, V, gs, dist_just)
        hard_astm = astm_jominy_curve(C, Mn, Si, Ni, Cr, Mo, Cu, gs, dist_just)
        exp_data = self.get_jominy_data()
        self.jominy_plot.plot_jominy(dist_just, hard_just, dist_just, hard_astm, exp_data, self.steel_combo.currentText(), hq, self.show_just_method)

    def save_composition(self):
        comp = self.get_composition()
        if comp is None: return
        file, _ = QFileDialog.getSaveFileName(self, "Save Composition", "", "JSON (*.json)")
        if file:
            with open(file, 'w') as f:
                json.dump({"composition": {"C":comp[0],"Mn":comp[1],"Si":comp[2],"Cr":comp[3],"Ni":comp[4],"Mo":comp[5],"Cu":comp[6],"V":comp[7],"W":comp[8]}, "grain_size":comp[9],"hardness_after_quenching":comp[10]}, f)
            QMessageBox.information(self, "Success", "Saved")

    def load_composition(self):
        file, _ = QFileDialog.getOpenFileName(self, "Load Composition", "", "JSON (*.json)")
        if file:
            with open(file, 'r') as f:
                data = json.load(f)
            comp = data.get("composition", {})
            self.c_input.setText(str(comp.get("C",0)))
            self.mn_input.setText(str(comp.get("Mn",0)))
            self.si_input.setText(str(comp.get("Si",0)))
            self.cr_input.setText(str(comp.get("Cr",0)))
            self.ni_input.setText(str(comp.get("Ni",0)))
            self.mo_input.setText(str(comp.get("Mo",0)))
            self.cu_input.setText(str(comp.get("Cu",0)))
            self.v_input.setText(str(comp.get("V",0)))
            self.w_input.setText(str(comp.get("W",0)))
            self.grain_size_input.setText(str(data.get("grain_size",7)))
            self.hardness_after_quenching_input.setText(str(data.get("hardness_after_quenching",50)))
            self.update_jominy_plot()

    def save_jominy_data(self):
        curves = self.get_jominy_data()
        if not curves:
            QMessageBox.warning(self,"No data")
            return
        file,_ = QFileDialog.getSaveFileName(self,"Save Jominy Data","","CSV (*.csv)")
        if file:
            d,h = curves[0]
            df = pd.DataFrame({"Distance":d, "Hardness":h})
            df.to_csv(file, index=False)

    def load_jominy_data(self):
        file,_ = QFileDialog.getOpenFileName(self,"Load Jominy Data","","CSV (*.csv)")
        if file:
            df = pd.read_csv(file)
            self.jominy_table.setRowCount(len(df))
            for i,row in df.iterrows():
                self.jominy_table.setItem(i,0, QTableWidgetItem(str(row.iloc[0])))
                self.jominy_table.setItem(i,1, QTableWidgetItem(str(row.iloc[1])))
            self.update_jominy_plot()

    def save_jominy_plot(self):
        file,_ = QFileDialog.getSaveFileName(self,"Save Plot","","PNG (*.png)")
        if file:
            self.jominy_plot.fig.savefig(file)

    def save_results(self):
        if not self.di_results:
            QMessageBox.warning(self,"No results")
            return
        file,_ = QFileDialog.getSaveFileName(self,"Save Results","","CSV (*.csv)")
        if file:
            with open(file,'w') as f:
                for k,v in self.di_results.items():
                    f.write(f"{k},{v}\n")
                if self.jominy_di:
                    f.write(f"Experimental actual,{self.jominy_di}\n")
                    f.write(f"Experimental ideal,{self.ideal_di}\n")

# =============================================================================
# 3. QUENCHING STUDIO – Steel Transformation Analysis (zwiększona czcionka, ikony)
# =============================================================================
steel_database_quenching = {
    "1.6580; 30CrNiMo8; 30H2N2M": {
        "C": 0.30, "Mn": 0.50, "Si": 0.30, "Cr": 2.00, "Ni": 2.00, "Mo": 0.30, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 120, "half_martensitic": 37,
        "description": "Heat treatable steel, used for shafts, axles, gears, high-strength bolts",
        "heat_treatment": "Quenching: 850-880°C in oil, Tempering: 550-650°C",
        "hardness": "280-320 HB",
        "hardness_after_quenching": "50-55 HRC"
    },
    "1.6582; 34CrNiMo6; 34HNM": {
        "C": 0.34, "Mn": 0.60, "Si": 0.30, "Cr": 1.50, "Ni": 1.50, "Mo": 0.25, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 100, "half_martensitic": 40,
        "description": "Heat treatable steel, used for high-strength machine components, crankshafts",
        "heat_treatment": "Quenching: 840-860°C in oil, Tempering: 580-650°C",
        "hardness": "300-340 HB",
        "hardness_after_quenching": "52-56 HRC"
    },
    "1.7035; 41Cr4; 40H": {
        "C": 0.41, "Mn": 0.70, "Si": 0.30, "Cr": 1.00, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 80, "half_martensitic": 43,
        "description": "Case-hardening and heat treatable steel, used for machine components, gears",
        "heat_treatment": "Quenching: 830-850°C in oil, Tempering: 550-650°C",
        "hardness": "200-250 HB (case-hardened), 250-300 HB (heat treated)",
        "hardness_after_quenching": "45-50 HRC"
    },
    "1.7225; 42CrMo4; 40HM": {
        "C": 0.42, "Mn": 0.75, "Si": 0.30, "Cr": 1.00, "Ni": 0.10, "Mo": 0.20, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 90, "half_martensitic": 43,
        "description": "Heat treatable steel, used for shafts, axles, fasteners, gears",
        "heat_treatment": "Quenching: 840-860°C in oil, Tempering: 550-650°C",
        "hardness": "280-320 HB",
        "hardness_after_quenching": "48-52 HRC"
    },
    "1.0503; C45; 45": {
        "C": 0.45, "Mn": 0.65, "Si": 0.25, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 30, "half_martensitic": 45,
        "description": "General purpose structural steel, used for machine components, shafts, axles",
        "heat_treatment": "Quenching: 820-850°C in water, Tempering: 550-600°C",
        "hardness": "170-220 HB",
        "hardness_after_quenching": "25-30 HRC"
    },
    "1.0535; C55; 55": {
        "C": 0.55, "Mn": 0.65, "Si": 0.25, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 35, "half_martensitic": 49,
        "description": "Higher strength structural steel, used for springs, machine components",
        "heat_treatment": "Quenching: 820-850°C in water, Tempering: 450-550°C",
        "hardness": "200-250 HB",
        "hardness_after_quenching": "30-35 HRC"
    },
    "1.0545; S355; 18G2A": {
        "C": 0.20, "Mn": 1.40, "Si": 0.40, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 25, "half_martensitic": 32,
        "description": "General purpose structural steel, used in construction and bridge building",
        "heat_treatment": "Normalizing: 890-950°C",
        "hardness": "140-180 HB",
        "hardness_after_quenching": "Not applicable (non-alloy steel)"
    },
    "1.0601; C60; 60": {
        "C": 0.60, "Mn": 0.75, "Si": 0.25, "Cr": 0.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 40, "half_martensitic": 51,
        "description": "Higher strength structural steel, used for springs, machine components",
        "heat_treatment": "Quenching: 810-840°C in water, Tempering: 400-500°C",
        "hardness": "220-270 HB",
        "hardness_after_quenching": "35-40 HRC"
    },
    "1.6587; 18CrNiMo7-6; 17HNM": {
        "C": 0.18, "Mn": 0.60, "Si": 0.30, "Cr": 1.70, "Ni": 1.50, "Mo": 0.30, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 110, "half_martensitic": 31,
        "description": "Case-hardening steel, used for gears, shafts, transmission components",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (surface after carburizing), 280-320 HB (core)",
        "hardness_after_quenching": "58-62 HRC (surface)"
    },
    "1.7131; 16MnCr5; 16HG": {
        "C": 0.16, "Mn": 1.10, "Si": 0.30, "Cr": 1.00, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 70, "half_martensitic": 30,
        "description": "Case-hardening steel, used for gears, shafts, machine components",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (surface after carburizing), 250-300 HB (core)",
        "hardness_after_quenching": "58-62 HRC (surface)"
    },
    "1.7147; 20MnCr5; 20HG": {
        "C": 0.20, "Mn": 1.10, "Si": 0.30, "Cr": 1.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 75, "half_martensitic": 32,
        "description": "Case-hardening steel, used for gears, shafts, machine components",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (surface after carburizing), 280-320 HB (core)",
        "hardness_after_quenching": "58-62 HRC (surface)"
    },
    "1.8509; 41CrAlMo7-10; 38HMJ": {
        "C": 0.41, "Mn": 0.60, "Si": 0.30, "Cr": 1.50, "Ni": 0.10, "Mo": 0.25, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 85, "half_martensitic": 43,
        "description": "Nitriding steel, used for components working at elevated temperatures",
        "heat_treatment": "Quenching: 930-950°C in oil, Tempering: 650-700°C",
        "hardness": "280-320 HB",
        "hardness_after_quenching": "45-50 HRC"
    },
}

def FahrenheitToCelsius(TF):
    return (TF - 32.)*5./9.

def FC(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return np.exp((1.0 + 6.31*C + 1.78*Mn + 0.31*Si + 1.12*Ni + 2.7*Cr + 4.06*Mo))

def PC(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return np.exp(-4.25 + 4.12*C + 4.36*Mn + 0.44*Si + 1.71*Ni + 3.33*Cr + 5.19*np.sqrt(Mo))

def BC(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return np.exp(-10.23 + 10.18*C + 0.85*Mn + 0.55*Ni + 0.9*Cr + 0.36*Mo)

def Ae1_Grange(**comp):
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    return FahrenheitToCelsius(1333 - 25*Mn + 40*Si - 26*Ni + 42*Cr)

def Ae3_Grange(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    return FahrenheitToCelsius(1570 - 323*C - 25*Mn + 80*Si - 32*Ni - 3*Cr)

def Ae1_Andrews(**comp):
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    W = comp.get('W', 0)
    As = comp.get('As', 0)
    return 723 - 16.9*Ni + 29.1*Si + 6.38*W - 10.7*Mn + 16.9*Cr + 290*As

def Ae3_Andrews(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    V = comp.get('V', 0)
    W = comp.get('W', 0)
    Cu = comp.get('Cu', 0)
    return 910 - 203*np.sqrt(C) + 44.7*Si - 15.2*Ni + 31.5*Mo + 104*V + 13.1*W - \
        30.0*Mn + 11.0*Cr + 20.0*Cu

def Bs_Li(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 637 - 58*C - 35*Mn - 15*Ni - 34*Cr - 41*Mo

def Bs_VanBohemen(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 839 - (86*Mn + 23*Si + 67*Cr + 33*Ni + 75*Mo) - 270*(1 - np.exp(-1.33*C))

def Ms_Andrews(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    Co = comp.get('Co', 0)
    return 539 - 423*C - 30.4*Mn - 17.7*Ni - 12.1*Cr - 7.5*Mo + 10*Co - 7.5*Si

def alpha_martensite_VanBohemen(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 1e-3*(27.2 - (0.14*Mn + 0.21*Si + 0.11*Cr + 0.08*Ni + 0.05*Mo) - 19.8*(1-np.exp(-1.56*C)))

def Ms_VanBohemen(**comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 565 - (31*Mn + 13*Si + 10*Cr + 18*Ni + 12*Mo) - 600*(1-np.exp(-0.96*C))

def Hv_martensite(phi700, **comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    return 127 + 949*C + 27*Si + 11*Mn + 8*Ni + 16*Cr + 21*np.log10(phi700*3600)

def Hv_bainite(phi700, **comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return -323 + 185*C + 330*Si + 153*Mn + 65*Ni + 144*Cr + 191*Mo + \
        (89 + 53*C - 55*Si - 22*Mn - 10*Ni - 20*Cr - 33*Mo)*np.log10(phi700*3600)

def Hv_ferrite_pearlite(phi700, **comp):
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    V = comp.get('V', 0)
    return 42 + 223*C + 53*Si + 30*Mn + 12.6*Ni + 7*Cr + 19*Mo + \
        (10 - 19*Si + 4*Ni + 8*Cr + 130*V)*np.log10(phi700*3600)

class Alloy:
    def __init__(self, gs, **w):
        self.gs = gs
        self.w = w
        self.C = w.get('C', 0)
        self.Mn = w.get('Mn', 0)
        self.Si = w.get('Si', 0)
        self.Ni = w.get('Ni', 0)
        self.Cr = w.get('Cr', 0)
        self.Mo = w.get('Mo', 0)
        self.Co = w.get('Co', 0)
        self.FC = FC(**w)
        self.PC = PC(**w)
        self.BC = BC(**w)
        self.Ae3 = Ae3_Andrews(**w)
        self.Ae1 = Ae1_Andrews(**w)
        self.Bs = Bs_Li(**w)
        self.Ms = Ms_VanBohemen(**w)
        self.alpha_martensite = alpha_martensite_VanBohemen(**w)
        self.Hv_martensite = lambda phi700: Hv_martensite(phi700, **w)
        self.Hv_bainite = lambda phi700: Hv_bainite(phi700, **w)
        self.Hv_ferrite_pearlite = lambda phi700: Hv_ferrite_pearlite(phi700, **w)

    def format_composition(self, vmin=0):
        fmt = []
        for k, v in self.w.items():
            if v > vmin:
                fmt.append('{:g}{:}'.format(v, k))
        fmt.insert(0, 'Fe')
        fmt = '-'.join(fmt) + ' (wt.%)'
        fmt += '\nASTM grain size {:}'.format(self.gs)
        return fmt

class SigmoidalFunction(object):
    tck = None
    tck_inv = None

    def __new__(cls, x):
        if cls is SigmoidalFunction:
            raise TypeError("Cannot instantiate abstract class SigmoidalFunction")
        for var in ['xmin', 'xmax', 'ymin', 'ymax', 'n']:
            if not hasattr(cls, var):
                raise NotImplementedError(
                    'Class {} does not have required attribute `{}`'.format(cls, var))
        return cls.val(x)

    @staticmethod
    def f(x):
        pass

    @classmethod
    def val(cls, x):
        if hasattr(x, '__iter__') and not isinstance(x, str):
            x = np.array(x)
            xmin = x[x > 0].min()
        else:
            xmin = x
        xmin = min(cls.xmin, xmin)
        if xmin < cls.xmin or cls.tck is None:
            cls.xmin = xmin
            cls.init_spline()
        return splev(x, cls.tck)

    @classmethod
    def inv(cls, y):
        if hasattr(y, '__iter__') and not isinstance(y, str):
            y = np.array(y)
            ymin, ymax = y.min(), y.max()
        else:
            ymin, ymax = y, y
        if cls.tck_inv is None:
            cls.init_spline()
        if ymin < cls.ymin or ymax > cls.ymax:
            print('Warning! y value out of range [{:g}:{:g}]. '
                  'Returned value is extrapolation'.format(cls.ymin, cls.ymax))
        return splev(y, cls.tck_inv)

    @classmethod
    def init_spline(cls):
        X = np.linspace(cls.xmin, cls.xmax, cls.n)
        Y = np.array([integrate.quad(cls.f, 0, x)[0] for x in X])
        cls.ymin = Y.min()
        cls.ymax = Y.max()
        cls.tck = splrep(X, Y)
        cls.tck_inv = splrep(Y, X)

class S(SigmoidalFunction):
    n = 999
    xmin = 0.001
    xmax = 0.999
    ymin = 0.02638507
    ymax = 2.02537893

    @staticmethod
    def f(x):
        return 1./(x**(0.4*(1. - x))*(1. - x)**(0.4*x))

class I(SigmoidalFunction):
    n = 999
    xmin = 0.001
    xmax = 0.999
    ymin = 0.29961765
    ymax = 4.05928646

    @staticmethod
    def f(x):
        return 1./(x**(2.*(1. - x)/3.)*(1. - x)**(2.*x/3.))

class PhaseTransformation(object):
    def __init__(self, alloy):
        self.alloy = alloy
        self.initialize()
        for var in ['comp_factor', 'Ts', 'Tf', 'Hv']:
            if not hasattr(self, var):
                raise NotImplementedError(
                    'Object {} does not have required attribute `{}`'.format(self, var))

    def __init_subclass__(cls):
        for var in ['Q', 'n1', 'n2']:
            if not hasattr(cls, var):
                raise NotImplementedError(
                    'Class {} does not have required attribute `{}`'.format(cls, var))

    def initialize(self):
        pass

    def get_transformation_factor(self, T):
        return self.comp_factor/(2**(self.n1*self.alloy.gs)*(self.Ts - T)**self.n2*np.exp(-self.Q/(R*(T + K))))

    def get_transformation_time(self, T, f):
        return S(f)*self.get_transformation_factor(T)

    def get_transformation_temperature(self, Tini, Tfin, cooling_rate, f, dT=1.0):
        dt = dT/np.array(cooling_rate)
        nt = len(dt) if hasattr(dt, '__len__') else 1
        T = np.arange(Tini, Tfin, -dT)
        nucleation_time = np.full((nt, len(T)), 0, dtype=float)
        filtr = T < self.Ts
        nucleation_time[:, filtr] = np.outer(dt, 1./self.get_transformation_factor(T[filtr]))
        nucleation_time = nucleation_time.cumsum(axis=1)
        Tt = np.full(nt, np.nan, dtype=float)
        Sf = S(f)
        for i, n_time in enumerate(nucleation_time):
            idx, = np.where(n_time >= Sf)
            if len(idx) > 0:
                Tt[i] = T[idx[0]]
        return float(Tt) if nt == 1 else Tt

    def get_transformed_fraction(self, t, T, n=1000):
        if len(t) > 3:
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            t2T = interp1d(t, T)
        dt = (max(t) - min(t))/(n - 1)
        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        nucleation_time = np.full(t.shape, 0, dtype=float)
        f = np.full(T.shape, 0, dtype=float)
        filtr = (T < self.Ts) & (T > self.Tf)
        if np.any(filtr):
            nucleation_time[filtr] = dt/self.get_transformation_factor(T[filtr])
            nucleation_time = nucleation_time.cumsum()
            if T[0] < self.Ts:
                nucleation_time += min(t)/self.get_transformation_factor(T[0])
            filtr = (nucleation_time >= S.ymin) & (nucleation_time <= S.ymax)
            if np.any(filtr):
                f[filtr] = S.inv(nucleation_time[filtr])
                f[nucleation_time < S.ymin] = 0
                f[nucleation_time > S.ymax] = 1
        return t, T, f

class Ferrite(PhaseTransformation):
    Q = 27500*4.184
    n1 = 0.41
    n2 = 3

    def initialize(self):
        self.comp_factor = self.alloy.FC
        self.Ts = self.alloy.Ae3
        self.Tf = self.alloy.Bs
        self.Hv = self.alloy.Hv_ferrite_pearlite

class Pearlite(PhaseTransformation):
    Q = 27500*4.184
    n1 = 0.32
    n2 = 3

    def initialize(self):
        self.comp_factor = self.alloy.PC
        self.Ts = self.alloy.Ae1
        self.Tf = self.alloy.Bs
        self.Hv = self.alloy.Hv_ferrite_pearlite

class Bainite(PhaseTransformation):
    Q = 27500*4.184
    n1 = 0.29
    n2 = 2

    def initialize(self):
        self.comp_factor = self.alloy.BC
        self.Ts = self.alloy.Bs
        self.Tf = self.alloy.Ms
        self.Hv = self.alloy.Hv_bainite

class Martensite:
    def __init__(self, alloy):
        self.alloy = alloy
        self.Ts = self.alloy.Ms
        self.Hv = self.alloy.Hv_martensite

    def get_transformed_fraction(self, t, T, n=1000):
        if len(t) > 3:
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            t2T = interp1d(t, T)
        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        f = np.full(T.shape, 0, dtype=float)
        filtr = T < self.alloy.Ms
        if np.any(filtr):
            f[filtr] = 1 - np.exp(-self.alloy.alpha_martensite*(self.alloy.Ms - T[filtr]))
        return t, T, f

class TransformationDiagrams:
    colors_dict = dict(ferrite='#1f77b4', pearlite='#ff7f0e', bainite='#2ca02c',
                       martensite='#d62728', austenite='#9467bd')
    columns_label_dict = dict(t='Time (s)', T='Temperature (°C)',
                              ferrite='Ferrite', pearlite='Pearlite', bainite='Bainite',
                              martensite='Martensite', austenite='Austenite')

    def __init__(self, alloy):
        self.alloy = alloy
        self.ferrite = Ferrite(self.alloy)
        self.pearlite = Pearlite(self.alloy)
        self.bainite = Bainite(self.alloy)
        self.martensite = Martensite(self.alloy)
        self.df_TTT = None
        self.df_CCT = None

    def get_transformed_fraction(self, t, T, n=1000):
        _, _, f_ferr = self.ferrite.get_transformed_fraction(t, T, n)
        _, _, f_pear = self.pearlite.get_transformed_fraction(t, T, n)
        _, _, f_bain = self.bainite.get_transformed_fraction(t, T, n)
        t, T, f_mart = self.martensite.get_transformed_fraction(t, T, n)
        f_ferr_inc = np.zeros(f_ferr.shape)
        f_pear_inc = np.zeros(f_pear.shape)
        f_bain_inc = np.zeros(f_bain.shape)
        f_mart_inc = np.zeros(f_mart.shape)
        f_ferr_inc[1:] = np.diff(f_ferr)
        f_pear_inc[1:] = np.diff(f_pear)
        f_bain_inc[1:] = np.diff(f_bain)
        f_mart_inc[1:] = np.diff(f_mart)
        f = pd.DataFrame({
            't': t.astype(float),
            'T': T.astype(float),
            'ferrite': np.zeros(len(t), dtype=float),
            'pearlite': np.zeros(len(t), dtype=float),
            'bainite': np.zeros(len(t), dtype=float),
            'martensite': np.zeros(len(t), dtype=float),
            'austenite': np.zeros(len(t), dtype=float)
        })
        f.loc[0, 'ferrite'] = float(f_ferr[0])
        f.loc[0, 'pearlite'] = float(f_pear[0])
        f.loc[0, 'bainite'] = float(f_bain[0])
        f.loc[0, 'martensite'] = float(f_mart[0])
        f.loc[0, 'austenite'] = float(1. - f_ferr[0] - f_pear[0] - f_bain[0] - f_mart[0])

        def f1(i, x, y, z, w):
            if f_ferr[i] < 1:
                return f.loc[i-1, 'ferrite'] + f_ferr_inc[i]*(1 - x - y - z - w)/(1 - f_ferr[i]) - x
            else:
                return f.loc[i-1, 'ferrite'] + f_ferr_inc[i]*(1 - y - z - w) - x

        def f2(i, x, y, z, w):
            if f_pear[i] < 1:
                return f.loc[i-1, 'pearlite'] + f_pear_inc[i]*(1 - x - y - z - w)/(1 - f_pear[i]) - y
            else:
                return f.loc[i-1, 'pearlite'] + f_pear_inc[i]*(1 - x - z - w) - y

        def f3(i, x, y, z, w):
            return f.loc[i-1, 'bainite'] + f_bain_inc[i]*(1 - x - y - w) - z

        def f4(i, x, y, z, w):
            return f.loc[i-1, 'martensite'] + f_mart_inc[i]*(1 - x - y - z) - w

        for i in range(1, len(f)):
            x0 = [f.loc[i-1, 'ferrite'], f.loc[i-1, 'pearlite'], f.loc[i-1, 'bainite'], f.loc[i-1, 'martensite']]
            res = root(lambda x: [f1(i, *x), f2(i, *x), f3(i, *x), f4(i, *x)], x0)
            f.loc[i, 'ferrite'] = float(res.x[0])
            f.loc[i, 'pearlite'] = float(res.x[1])
            f.loc[i, 'bainite'] = float(res.x[2])
            f.loc[i, 'martensite'] = float(res.x[3])
            f.loc[i, 'austenite'] = float(1.0 - sum(res.x))

        phi700 = None
        try:
            T2t = interp1d(T, t)
            phi700 = 2./(T2t(699.) - T2t(701.))
            if phi700 == 0:
                phi700 = None
        except ValueError:
            pass
        if phi700 is not None:
            f['Hv'] = f['martensite']*self.martensite.Hv(phi700) + f['bainite']*self.bainite.Hv(phi700) + \
                (f['ferrite'] + f['pearlite'])*self.ferrite.Hv(phi700)
        else:
            f['Hv'] = np.nan
        return f.round(12)

    def draw_thermal_cycle(self, ax, t, T, n=100, **kwargs):
        if len(t) > 3:
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            t2T = interp1d(t, T)
        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        kw = dict(color='k', ls='--')
        kw.update(kwargs)
        return ax.plot(t, T, **kw)

    def TTT(self, fs=1e-2, ff=.99, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()
        T = np.arange(self.alloy.Bs, self.alloy.Ae3)
        ts = self.ferrite.get_transformation_time(T, fs)
        tf = self.ferrite.get_transformation_time(T, ff)
        ax.plot(ts, T, color=self.colors_dict['ferrite'], label='Ferrite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['ferrite'], ls='--', label='Ferrite {:.0f}%'.format(100*ff), **kwargs)
        T = np.arange(self.alloy.Bs, self.alloy.Ae1)
        ts = self.pearlite.get_transformation_time(T, fs)
        tf = self.pearlite.get_transformation_time(T, ff)
        ax.plot(ts, T, color=self.colors_dict['pearlite'], label='Pearlite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['pearlite'], ls='--', label='Pearlite {:.0f}%'.format(100*ff), **kwargs)
        T = np.arange(self.alloy.Ms, self.alloy.Bs)
        ts = self.bainite.get_transformation_time(T, fs)
        tf = self.bainite.get_transformation_time(T, ff)
        ax.plot(ts, T, color=self.colors_dict['bainite'], label='Bainite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['bainite'], ls='--', label='Bainite {:.0f}%'.format(100*ff), **kwargs)
        ax.axhline(self.alloy.Ae3, xmax=.1, color=self.colors_dict['ferrite'], ls=':')
        ax.axhline(self.alloy.Ae1, xmax=.1, color=self.colors_dict['pearlite'], ls=':')
        ax.axhline(self.alloy.Bs, color=self.colors_dict['bainite'], ls=':')
        ax.axhline(self.alloy.Ms, color=self.colors_dict['martensite'])
        ax.set_xscale('log')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Temperature (°C)')
        ax.set_title(self.alloy.format_composition())
        xmin = ax.get_xlim()[0]
        ax.text(xmin*1.5, self.alloy.Ae3, 'Ae3', color=self.colors_dict['ferrite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ae1, 'Ae1', color=self.colors_dict['pearlite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Bs, 'Bs', color=self.colors_dict['bainite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ms, 'Ms', color=self.colors_dict['martensite'], ha='left', va='bottom')
        ax.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, -.15))
        fig.subplots_adjust(bottom=.2)
        return ax

    def CCT(self, Tini=900, fs=1e-2, ff=.99, phi_min=1e-4, phi_max=1e4, phi_steps=420, ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()
        cooling_rates = 10**np.linspace(np.log10(phi_min), np.log10(phi_max), phi_steps)
        draw_cooling = kwargs.get('draw_cooling', True)
        Ts = self.ferrite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, fs)
        Tf = self.ferrite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, ff)
        ax.plot(Ts/cooling_rates, Ts, color=self.colors_dict['ferrite'], label='Ferrite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['ferrite'], ls='--', label='Ferrite {:.0f}%'.format(100*ff), **kwargs)
        Ts = self.pearlite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, fs)
        Tf = self.pearlite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, ff)
        ax.plot(Ts/cooling_rates, Ts, color=self.colors_dict['pearlite'], label='Pearlite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['pearlite'], ls='--', label='Pearlite {:.0f}%'.format(100*ff), **kwargs)
        Ts = self.bainite.get_transformation_temperature(Tini, self.alloy.Ms, cooling_rates, fs)
        Tf = self.bainite.get_transformation_temperature(Tini, self.alloy.Ms, cooling_rates, ff)
        ax.plot(Ts/cooling_rates, Ts, color=self.colors_dict['bainite'], label='Bainite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['bainite'], ls='--', label='Bainite {:.0f}%'.format(100*ff), **kwargs)
        ax.axhline(self.alloy.Ae3, xmax=.1, color=self.colors_dict['ferrite'], ls=':')
        ax.axhline(self.alloy.Ae1, xmax=.1, color=self.colors_dict['pearlite'], ls=':')
        ax.axhline(self.alloy.Bs, color=self.colors_dict['bainite'], ls=':')
        ax.axhline(self.alloy.Ms, color=self.colors_dict['martensite'])
        if draw_cooling:
            for cooling_rate in cooling_rates[::10]:
                T = np.linspace(Tini, 25, 100)
                t = (Tini - T)/cooling_rate
                kw = dict(lw=.5)
                kw.update(kwargs)
                ax.plot(t, T, 'k:', **kw)
        ax.set_xscale('log')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Temperature (°C)')
        ax.set_title(self.alloy.format_composition())
        xmin = ax.get_xlim()[0]
        ax.text(xmin*1.5, self.alloy.Ae3, 'Ae3', color=self.colors_dict['ferrite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ae1, 'Ae1', color=self.colors_dict['pearlite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Bs, 'Bs', color=self.colors_dict['bainite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ms, 'Ms', color=self.colors_dict['martensite'], ha='left', va='bottom')
        ax.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, -.15))
        fig.subplots_adjust(bottom=.2)
        return ax

    def plot_phase_fraction(self, t, T, n=1000, xaxis='t', ax=None, **kwargs):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        if len(t) > 3:
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            t2T = interp1d(t, T)
        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        f = self.get_transformed_fraction(t, T, n)
        if f['ferrite'].max() > 0:
            ax.plot(f[xaxis], f['ferrite'], color=self.colors_dict['ferrite'], label='Ferrite')
        if f['pearlite'].max() > 0:
            ax.plot(f[xaxis], f['pearlite'], color=self.colors_dict['pearlite'], label='Pearlite')
        if f['bainite'].max() > 0:
            ax.plot(f[xaxis], f['bainite'], color=self.colors_dict['bainite'], label='Bainite')
        if f['martensite'].max() > 0:
            ax.plot(f[xaxis], f['martensite'], color=self.colors_dict['martensite'], label='Martensite')
        if f['austenite'].max() > 0:
            ax.plot(f[xaxis], f['austenite'], color=self.colors_dict['austenite'], label='Austenite')
        if not np.isnan(f.iloc[-1]['Hv']):
            T_ref = 25
            try:
                Hv_ref = interp1d(f['T'], f['Hv'])(T_ref)
            except ValueError:
                T_ref, Hv_ref = f.iloc[-1]['T'], f.iloc[-1]['Hv']
            ax.text(.95, .95, 'Hardness for phase fractions at {:.1f} °C: {:.0f} HV'.format(T_ref, Hv_ref),
                    transform=ax.transAxes, ha='right', va='top',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        ax.set_xlabel(self.columns_label_dict[xaxis])
        ax.set_ylabel('Phase Fraction')
        ax.legend()
        return ax

class SteelAnalysisApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Quenching studio 1.2")
        self.setGeometry(100, 50, 1800, 1000)
        # Zwiększenie czcionki
        font = QFont("Segoe UI", 10)
        self.setFont(font)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        self.input_tab = QWidget()
        input_layout = QHBoxLayout(self.input_tab)
        col1_layout = QVBoxLayout()
        steel_group = QGroupBox("Steel Selection")
        steel_layout = QVBoxLayout(steel_group)
        self.steel_combo = QComboBox()
        self.steel_combo.addItems(steel_database_quenching.keys())
        self.steel_combo.currentTextChanged.connect(self.load_steel_data)
        steel_layout.addWidget(QLabel("Select steel:"))
        steel_layout.addWidget(self.steel_combo)
        self.steel_info = QTextEdit()
        self.steel_info.setReadOnly(True)
        self.steel_info.setMaximumHeight(200)
        steel_layout.addWidget(QLabel("Steel information:"))
        steel_layout.addWidget(self.steel_info)
        col1_layout.addWidget(steel_group)
        col1_layout.addStretch()
        col2_layout = QVBoxLayout()
        composition_group = QGroupBox("Alloy Composition (wt.%)")
        composition_layout = QGridLayout(composition_group)
        self.composition_inputs = {}
        elements = ['C', 'Mn', 'Si', 'Ni', 'Cr', 'Mo', 'V', 'Co', 'Cu', 'W']
        defaults = [0.37, 0.77, 0.15, 0.04, 0.98, 0.21, 0.0, 0.0, 0.0, 0.0]
        for i, element in enumerate(elements):
            composition_layout.addWidget(QLabel(element), i, 0)
            spinbox = QDoubleSpinBox()
            spinbox.setRange(0, 100)
            spinbox.setValue(defaults[i])
            spinbox.setDecimals(3)
            spinbox.setSingleStep(0.01)
            self.composition_inputs[element] = spinbox
            composition_layout.addWidget(spinbox, i, 1)
        save_load_layout = QHBoxLayout()
        self.save_composition_btn = QPushButton(" Save composition")
        self.save_composition_btn.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_composition_btn.clicked.connect(self.save_composition)
        self.load_composition_btn = QPushButton(" Load composition")
        self.load_composition_btn.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.load_composition_btn.clicked.connect(self.load_composition)
        save_load_layout.addWidget(self.save_composition_btn)
        save_load_layout.addWidget(self.load_composition_btn)
        composition_layout.addLayout(save_load_layout, len(elements), 0, 1, 2)
        composition_layout.addWidget(QLabel("Grain Size"), len(elements)+1, 0)
        self.gs_input = QDoubleSpinBox()
        self.gs_input.setRange(1, 12)
        self.gs_input.setValue(7)
        self.gs_input.setDecimals(1)
        composition_layout.addWidget(self.gs_input, len(elements)+1, 1)
        col2_layout.addWidget(composition_group)
        col2_layout.addStretch()
        col3_layout = QVBoxLayout()
        temp_group = QGroupBox("Temperature Parameters (°C)")
        temp_layout = QGridLayout(temp_group)
        temp_layout.addWidget(QLabel("Initial Temperature (Tini)"), 0, 0)
        self.tini_input = QDoubleSpinBox()
        self.tini_input.setRange(0, 2000)
        self.tini_input.setValue(900)
        temp_layout.addWidget(self.tini_input, 0, 1)
        temp_layout.addWidget(QLabel("Final Temperature (Tfin)"), 1, 0)
        self.tfin_input = QDoubleSpinBox()
        self.tfin_input.setRange(0, 1000)
        self.tfin_input.setValue(25)
        temp_layout.addWidget(self.tfin_input, 1, 1)
        temp_layout.addWidget(QLabel("Single Cooling Rate (°C/s)"), 2, 0)
        self.single_cooling_input = QDoubleSpinBox()
        self.single_cooling_input.setRange(0.001, 10000)
        self.single_cooling_input.setValue(10)
        self.single_cooling_input.setDecimals(3)
        self.single_cooling_input.setSingleStep(1)
        temp_layout.addWidget(self.single_cooling_input, 2, 1)
        col3_layout.addWidget(temp_group)
        cooling_group = QGroupBox("Cooling Rates for Analysis (°C/s)")
        cooling_layout = QVBoxLayout(cooling_group)
        self.cooling_rates_edit = QTextEdit()
        self.cooling_rates_edit.setPlainText(
            "1000\n300\n100\n30\n10\n3\n1\n0.3\n0.1\n0.03\n0.01\n0.003\n0.001"
        )
        self.cooling_rates_edit.setMaximumHeight(150)
        cooling_layout.addWidget(self.cooling_rates_edit)
        col3_layout.addWidget(cooling_group)
        
        self.calculate_single_btn = QPushButton(" Calculate Single Cooling Rate")
        self.calculate_single_btn.setIcon(self.style().standardIcon(QStyle.SP_ComputerIcon))
        self.calculate_single_btn.clicked.connect(self.calculate_single_cooling)
        self.calculate_hardness_btn = QPushButton(" Calculate Hardness vs Cooling Rate")
        self.calculate_hardness_btn.setIcon(self.style().standardIcon(QStyle.SP_FileDialogDetailedView))
        self.calculate_hardness_btn.clicked.connect(self.calculate_hardness)
        self.calculate_ttt_btn = QPushButton(" Calculate TTT Diagram")
        self.calculate_ttt_btn.setIcon(self.style().standardIcon(QStyle.SP_FileDialogContentsView))
        self.calculate_ttt_btn.clicked.connect(self.calculate_ttt)
        self.calculate_cct_btn = QPushButton(" Calculate CCT Diagram")
        self.calculate_cct_btn.setIcon(self.style().standardIcon(QStyle.SP_FileDialogInfoView))
        self.calculate_cct_btn.clicked.connect(self.calculate_cct)
        col3_layout.addWidget(self.calculate_single_btn)
        col3_layout.addWidget(self.calculate_hardness_btn)
        col3_layout.addWidget(self.calculate_ttt_btn)
        col3_layout.addWidget(self.calculate_cct_btn)
        col3_layout.addStretch()
        input_layout.addLayout(col1_layout)
        input_layout.addLayout(col2_layout)
        input_layout.addLayout(col3_layout)
        self.tab_widget.addTab(self.input_tab, "Input Parameters")
        self.single_tab = QWidget()
        single_layout = QVBoxLayout(self.single_tab)
        self.single_fig = Figure(figsize=(10, 8))
        self.single_canvas = FigureCanvas(self.single_fig)
        self.single_toolbar = NavigationToolbar(self.single_canvas, self)
        single_layout.addWidget(self.single_toolbar)
        single_layout.addWidget(self.single_canvas)
        self.single_table_tab = QWidget()
        single_table_layout = QVBoxLayout(self.single_table_tab)
        self.single_table = QTableWidget()
        self.single_table.setColumnCount(2)
        self.single_table.setHorizontalHeaderLabels(["Property", "Value"])
        self.single_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        single_table_layout.addWidget(self.single_table)
        self.save_single_table_btn = QPushButton(" Save results to CSV")
        self.save_single_table_btn.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_single_table_btn.clicked.connect(self.save_single_table)
        single_table_layout.addWidget(self.save_single_table_btn)
        self.hardness_tab = QWidget()
        hardness_layout = QVBoxLayout(self.hardness_tab)
        self.fig1 = Figure(figsize=(10, 8))
        self.canvas1 = FigureCanvas(self.fig1)
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        hardness_layout.addWidget(self.toolbar1)
        hardness_layout.addWidget(self.canvas1)
        self.hardness_table_tab = QWidget()
        hardness_table_layout = QVBoxLayout(self.hardness_table_tab)
        self.hardness_table = QTableWidget()
        self.hardness_table.setColumnCount(6)
        self.hardness_table.setHorizontalHeaderLabels(["Cooling Rate (°C/s)", "Ferrite", "Pearlite", "Bainite", "Martensite", "Hardness (HV)"])
        self.hardness_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        hardness_table_layout.addWidget(self.hardness_table)
        self.save_hardness_table_btn = QPushButton(" Save results to CSV")
        self.save_hardness_table_btn.setIcon(self.style().standardIcon(QStyle.SP_DialogSaveButton))
        self.save_hardness_table_btn.clicked.connect(self.save_hardness_table)
        hardness_table_layout.addWidget(self.save_hardness_table_btn)
        self.ttt_tab = QWidget()
        ttt_layout = QVBoxLayout(self.ttt_tab)
        self.fig2 = Figure(figsize=(10, 8))
        self.canvas2 = FigureCanvas(self.fig2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        ttt_layout.addWidget(self.toolbar2)
        ttt_layout.addWidget(self.canvas2)
        self.cct_tab = QWidget()
        cct_layout = QVBoxLayout(self.cct_tab)
        self.fig3 = Figure(figsize=(10, 8))
        self.canvas3 = FigureCanvas(self.fig3)
        self.toolbar3 = NavigationToolbar(self.canvas3, self)
        cct_layout.addWidget(self.toolbar3)
        cct_layout.addWidget(self.canvas3)
        self.instruction_tab = QWidget()
        instruction_layout = QVBoxLayout(self.instruction_tab)
        instruction_text = QTextEdit()
        instruction_text.setReadOnly(True)
        instruction_text.setPlainText("""
This application is based on GNU code from the repository:
        https://github.com/arthursn/transformation-diagrams/tree/master

The program has educational and illustrative character for students
The author is not responsible for errors and other effects of the program
(c) Marek Góral 2025 m_goral@interia.pl
Program created as part of the Heat Treatment and Thermochemical Treatment course
within the Regional Excellence Initiative program at Rzeszow University of Technology
Free program for educational and didactic use

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

STEEL TRANSFORMATION ANALYSIS SOFTWARE USER MANUAL

1. SELECT STEEL
   - Select steel from the dropdown list in the "Input Parameters" tab
   - The program will automatically load chemical composition and steel parameters

2. ADJUST PARAMETERS (OPTIONAL)
   - You can manually change chemical composition in input fields
   - Adjust ASTM grain size
   - Set initial and final temperature
   - Enter desired cooling rate or list of cooling rates

3. CALCULATIONS
   - Click appropriate calculation button:
     * "Calculate Single Cooling Rate" - for one cooling rate
     * "Calculate Hardness vs Cooling Rate" - for multiple cooling rates
     * "Calculate TTT Diagram" - isothermal transformation diagram
     * "Calculate CCT Diagram" - continuous cooling transformation diagram

4. RESULTS
   - Results are displayed in appropriate tabs:
     * "Single Cooling Rate" - phase fraction charts
     * "Single Cooling Rate - Table" - results table for one cooling rate
     * "Hardness vs Cooling Rate" - dependency charts
     * "Hardness vs Cooling Rate - Table" - results table for multiple cooling rates
     * "TTT Diagram" - isothermal transformation diagram
     * "CCT Diagram" - continuous cooling transformation diagram

5. SAVING RESULTS
   - In table tabs use "Save results to CSV" button
   - You can save chemical composition with "Save composition" button

PROGRAM DESCRIPTION:

The program is used for analysis of phase transformations in steels during heat treatment.
Based on chemical composition of steel it calculates:
- Critical temperatures (Ae1, Ae3, Bs, Ms)
- Kinetics of phase transformations (ferritic, pearlitic, bainitic, martensitic)
- Final hardness as function of cooling rate
- TTT diagrams (Time-Temperature-Transformation)
- CCT diagrams (Continuous-Cooling-Transformation)

The program uses empirical models based on scientific research
to predict steel behavior during various heat treatment processes.
        """)
        instruction_layout.addWidget(instruction_text)
        self.tab_widget.addTab(self.single_tab, "Single Cooling Rate")
        self.tab_widget.addTab(self.single_table_tab, "Single Cooling Rate - Table")
        self.tab_widget.addTab(self.hardness_tab, "Hardness vs Cooling Rate")
        self.tab_widget.addTab(self.hardness_table_tab, "Hardness vs Cooling Rate - Table")
        self.tab_widget.addTab(self.ttt_tab, "TTT Diagram")
        self.tab_widget.addTab(self.cct_tab, "CCT Diagram")
        self.tab_widget.addTab(self.instruction_tab, "Instructions")
        self.alloy = None
        self.diagrams = None
        self.single_results = None
        self.hardness_results = None
        self.load_steel_data()

    def load_steel_data(self):
        steel_name = self.steel_combo.currentText()
        steel_data = steel_database_quenching[steel_name]
        for element, value in steel_data.items():
            if element in self.composition_inputs:
                self.composition_inputs[element].setValue(value)
        self.gs_input.setValue(steel_data.get("grain_size", 7))
        info_text = f"Name: {steel_name}\n\n"
        info_text += f"Description: {steel_data.get('description', 'No description')}\n\n"
        info_text += f"Heat Treatment: {steel_data.get('heat_treatment', 'No heat treatment information')}\n\n"
        info_text += f"Typical Hardness: {steel_data.get('hardness', 'No hardness information')}\n"
        info_text += f"Hardness after Quenching: {steel_data.get('hardness_after_quenching', 'No hardness after quenching information')}\n"
        info_text += f"Typical DI: {steel_data.get('typical_DI', 'N/A')}\n"
        info_text += f"Half Martensitic: {steel_data.get('half_martensitic', 'N/A')}°C"
        self.steel_info.setPlainText(info_text)

    def get_composition(self):
        composition = {}
        for element, spinbox in self.composition_inputs.items():
            composition[element] = spinbox.value()
        return composition

    def save_composition(self):
        try:
            options = QFileDialog.Options()
            file_name, _ = QFileDialog.getSaveFileName(
                self, "Save Composition", "", "JSON Files (*.json)", options=options)
            if file_name:
                if not file_name.endswith('.json'):
                    file_name += '.json'
                composition = self.get_composition()
                data = {'composition': composition, 'grain_size': self.gs_input.value()}
                with open(file_name, 'w') as f:
                    json.dump(data, f, indent=4)
                QMessageBox.information(self, "Success", "Composition saved successfully!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save composition: {str(e)}")

    def load_composition(self):
        try:
            options = QFileDialog.Options()
            file_name, _ = QFileDialog.getOpenFileName(
                self, "Load Composition", "", "JSON Files (*.json)", options=options)
            if file_name:
                with open(file_name, 'r') as f:
                    data = json.load(f)
                composition = data.get('composition', {})
                grain_size = data.get('grain_size', 7)
                for element, value in composition.items():
                    if element in self.composition_inputs:
                        self.composition_inputs[element].setValue(value)
                self.gs_input.setValue(grain_size)
                self.steel_combo.setCurrentText("Custom")
                QMessageBox.information(self, "Success", "Composition loaded successfully!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load composition: {str(e)}")

    def calculate_single_cooling(self):
        try:
            composition = self.get_composition()
            gs = self.gs_input.value()
            Tini = self.tini_input.value()
            Tfin = self.tfin_input.value()
            cooling_rate = self.single_cooling_input.value()
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            total_time = (Tini - Tfin) / cooling_rate
            f = self.diagrams.get_transformed_fraction([0, total_time], [Tini, Tfin])
            f_fin = f.iloc[-1]
            self.single_results = {
                "Cooling Rate": cooling_rate,
                "Ferrite": f_fin['ferrite'],
                "Pearlite": f_fin['pearlite'],
                "Bainite": f_fin['bainite'],
                "Martensite": f_fin['martensite'],
                "Hardness": f_fin['Hv']
            }
            self.update_single_cooling_plot(f, cooling_rate)
            self.update_single_results_table(f_fin, cooling_rate)
            self.tab_widget.setCurrentIndex(1)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")

    def update_single_cooling_plot(self, f, cooling_rate):
        self.single_fig.clear()
        ax1 = self.single_fig.add_subplot(121)
        ax2 = self.single_fig.add_subplot(122)
        phases = ['ferrite', 'pearlite', 'bainite', 'martensite', 'austenite']
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        for phase, color in zip(phases, colors):
            if f[phase].max() > 0:
                ax1.plot(f['t'], f[phase], color=color, label=phase.capitalize())
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Phase Fraction')
        ax1.set_title(f'Phase Fractions vs Time (Cooling Rate: {cooling_rate} °C/s)')
        ax1.legend()
        ax1.grid(True)
        phase_fractions = [f.iloc[-1][phase] for phase in phases]
        labels = [f'{phase.capitalize()}: {frac:.1%}' for phase, frac in zip(phases, phase_fractions) if frac > 0]
        sizes = [frac for frac in phase_fractions if frac > 0]
        colors_pie = [color for phase, color in zip(phases, colors) if f.iloc[-1][phase] > 0]
        if sizes:
            ax2.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%', startangle=90)
            ax2.axis('equal')
            ax2.set_title('Final Phase Fractions')
        self.single_fig.suptitle(self.alloy.format_composition())
        self.single_canvas.draw()

    def update_single_results_table(self, f_fin, cooling_rate):
        self.single_table.setRowCount(6)
        self.single_table.setColumnCount(2)
        self.single_table.setHorizontalHeaderLabels(["Property", "Value"])
        data = [
            ["Cooling Rate", f"{cooling_rate} °C/s"],
            ["Ferrite", f"{f_fin['ferrite']:.3f}"],
            ["Pearlite", f"{f_fin['pearlite']:.3f}"],
            ["Bainite", f"{f_fin['bainite']:.3f}"],
            ["Martensite", f"{f_fin['martensite']:.3f}"],
            ["Hardness", f"{f_fin['Hv']:.0f} HV"]
        ]
        for row, (prop, value) in enumerate(data):
            self.single_table.setItem(row, 0, QTableWidgetItem(prop))
            self.single_table.setItem(row, 1, QTableWidgetItem(value))

    def save_single_table(self):
        try:
            options = QFileDialog.Options()
            file_name, _ = QFileDialog.getSaveFileName(
                self, "Save Results", "", "CSV Files (*.csv)", options=options)
            if file_name:
                if not file_name.endswith('.csv'):
                    file_name += '.csv'
                data = []
                for row in range(self.single_table.rowCount()):
                    prop = self.single_table.item(row, 0).text()
                    value = self.single_table.item(row, 1).text()
                    data.append([prop, value])
                df = pd.DataFrame(data, columns=["Property", "Value"])
                df.to_csv(file_name, index=False, encoding='utf-8')
                QMessageBox.information(self, "Success", "Results saved successfully!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save results: {str(e)}")

    def calculate_hardness(self):
        try:
            composition = self.get_composition()
            gs = self.gs_input.value()
            Tini = self.tini_input.value()
            Tfin = self.tfin_input.value()
            cooling_text = self.cooling_rates_edit.toPlainText()
            cooling_rates = []
            for line in cooling_text.split('\n'):
                try:
                    cooling_rates.append(float(line.strip()))
                except ValueError:
                    pass
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            f_ferr, f_pear, f_bain, f_mart, Hv = [], [], [], [], []
            for phi in cooling_rates:
                total_time = (Tini - Tfin) / phi
                f = self.diagrams.get_transformed_fraction([0, total_time], [Tini, Tfin])
                f_fin = f.iloc[-1]
                f_ferr.append(f_fin['ferrite'])
                f_pear.append(f_fin['pearlite'])
                f_bain.append(f_fin['bainite'])
                f_mart.append(f_fin['martensite'])
                Hv.append(f_fin['Hv'])
            self.hardness_results = pd.DataFrame({
                "Cooling Rate": cooling_rates,
                "Ferrite": f_ferr,
                "Pearlite": f_pear,
                "Bainite": f_bain,
                "Martensite": f_mart,
                "Hardness": Hv
            })
            self.update_hardness_plots(cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv)
            self.update_results_table(cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv)
            self.tab_widget.setCurrentIndex(3)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")

    def update_hardness_plots(self, cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv):
        self.fig1.clear()
        ax1 = self.fig1.add_subplot(121)
        ax2 = self.fig1.add_subplot(122)
        ax1.plot(cooling_rates, f_ferr, label='Ferrite')
        ax1.plot(cooling_rates, f_pear, label='Pearlite')
        ax1.plot(cooling_rates, f_bain, label='Bainite')
        ax1.plot(cooling_rates, f_mart, label='Martensite')
        ax1.set_xlabel('Cooling Rate (°C/s)')
        ax1.set_ylabel('Phase Fraction')
        ax1.set_xscale('log')
        ax1.set_title('Phase Fractions vs Cooling Rate')
        ax1.legend()
        ax1.grid(True)
        ax2.plot(cooling_rates, Hv, 'b-')
        ax2.set_xlabel('Cooling Rate (°C/s)')
        ax2.set_ylabel('Vickers Hardness (HV)')
        ax2.set_xscale('log')
        ax2.set_title('Hardness vs Cooling Rate')
        ax2.grid(True)
        self.fig1.suptitle(self.alloy.format_composition())
        self.canvas1.draw()

    def update_results_table(self, cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv):
        self.hardness_table.setRowCount(len(cooling_rates))
        self.hardness_table.setColumnCount(6)
        self.hardness_table.setHorizontalHeaderLabels(["Cooling Rate (°C/s)", "Ferrite", "Pearlite", "Bainite", "Martensite", "Hardness (HV)"])
        for i, (cr, ferr, pear, bain, mart, hv) in enumerate(zip(cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv)):
            self.hardness_table.setItem(i, 0, QTableWidgetItem(f"{cr:.3f}"))
            self.hardness_table.setItem(i, 1, QTableWidgetItem(f"{ferr:.3f}"))
            self.hardness_table.setItem(i, 2, QTableWidgetItem(f"{pear:.3f}"))
            self.hardness_table.setItem(i, 3, QTableWidgetItem(f"{bain:.3f}"))
            self.hardness_table.setItem(i, 4, QTableWidgetItem(f"{mart:.3f}"))
            self.hardness_table.setItem(i, 5, QTableWidgetItem(f"{hv:.0f}"))

    def save_hardness_table(self):
        try:
            options = QFileDialog.Options()
            file_name, _ = QFileDialog.getSaveFileName(
                self, "Save Results", "", "CSV Files (*.csv)", options=options)
            if file_name:
                if not file_name.endswith('.csv'):
                    file_name += '.csv'
                data = []
                for row in range(self.hardness_table.rowCount()):
                    row_data = []
                    for col in range(self.hardness_table.columnCount()):
                        item = self.hardness_table.item(row, col)
                        if item is not None:
                            row_data.append(item.text())
                        else:
                            row_data.append("")
                    data.append(row_data)
                df = pd.DataFrame(data, columns=["Cooling Rate (°C/s)", "Ferrite", "Pearlite", "Bainite", "Martensite", "Hardness (HV)"])
                df.to_csv(file_name, index=False, encoding='utf-8')
                QMessageBox.information(self, "Success", "Results saved successfully!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save results: {str(e)}")

    def calculate_ttt(self):
        try:
            composition = self.get_composition()
            gs = self.gs_input.value()
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            self.update_ttt_diagram()
            self.tab_widget.setCurrentIndex(5)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")

    def update_ttt_diagram(self):
        self.fig2.clear()
        ax = self.fig2.add_subplot(111)
        self.diagrams.TTT(ax=ax)
        self.canvas2.draw()

    def calculate_cct(self):
        try:
            composition = self.get_composition()
            gs = self.gs_input.value()
            Tini = self.tini_input.value()
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            self.update_cct_diagram(Tini)
            self.tab_widget.setCurrentIndex(6)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")

    def update_cct_diagram(self, Tini):
        self.fig3.clear()
        ax = self.fig3.add_subplot(111)
        self.diagrams.CCT(Tini=Tini, ax=ax)
        self.canvas3.draw()

# =============================================================================
# MAIN LAUNCHER (z nową kolejnością, przyciskami Close Window i About)
# =============================================================================
class SteelSuiteLauncher(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Heat treatment edu pack 1.2")
        self.setGeometry(300, 300, 550, 450)
        # Zwiększenie czcionki w launcherze
        font = QFont("Segoe UI", 11)
        self.setFont(font)
        
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setSpacing(15)
        
        title = QLabel("🔥 Heat treatment edu pack 1.2")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("font-size: 20pt; font-weight: bold; margin: 15px; color: #2c3e50;")
        layout.addWidget(title)
        
        # Kolejność: Quenching Studio, Jomina Analyser, Temperdata
        self.btn_quench = QPushButton("💧 Quenching Studio (Phase Transformations)")
        self.btn_quench.setIcon(self.style().standardIcon(QStyle.SP_ComputerIcon))
        self.btn_quench.clicked.connect(self.open_quenching)
        
        self.btn_jom = QPushButton("📊 Jomina Analyser (Hardenability)")
        self.btn_jom.setIcon(self.style().standardIcon(QStyle.SP_FileDialogDetailedView))
        self.btn_jom.clicked.connect(self.open_jomina)
        
        self.btn_temp = QPushButton("🌡️ Temperdata (Tempering Analysis)")
        self.btn_temp.setIcon(self.style().standardIcon(QStyle.SP_FileDialogContentsView))
        self.btn_temp.clicked.connect(self.open_temperdata)
        
        # Dodatkowe przyciski
        self.btn_close_window = QPushButton("❌ Close Window")
        self.btn_close_window.setIcon(self.style().standardIcon(QStyle.SP_DialogCloseButton))
        self.btn_close_window.clicked.connect(self.close_active_window)
        
        self.btn_about = QPushButton("ℹ️ About")
        self.btn_about.setIcon(self.style().standardIcon(QStyle.SP_MessageBoxInformation))
        self.btn_about.clicked.connect(self.show_about)
        
        # Style dla przycisków
        button_style = """
            QPushButton {
                font-size: 12pt;
                padding: 8px;
                text-align: left;
                background-color: #ecf0f1;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #bdc3c7;
            }
        """
        for btn in [self.btn_quench, self.btn_jom, self.btn_temp, self.btn_close_window, self.btn_about]:
            btn.setStyleSheet(button_style)
            layout.addWidget(btn)
        
        layout.addStretch()
        self.status_label = QLabel("Ready")
        layout.addWidget(self.status_label)
        
        self.windows = []  # przechowuje referencje do otwartych okien

    def open_temperdata(self):
        self.win = TemperingCurveApp()
        self.win.show()
        self.windows.append(self.win)

    def open_jomina(self):
        self.win = SteelCalculatorApp()
        self.win.show()
        self.windows.append(self.win)

    def open_quenching(self):
        self.win = SteelAnalysisApp()
        self.win.show()
        self.windows.append(self.win)

    def close_active_window(self):
        if self.windows:
            w = self.windows.pop()
            w.close()
            self.status_label.setText("Closed a window.")
        else:
            QMessageBox.information(self, "No window", "There is no open child window.")

    def show_about(self):
        about_text = """
        <h2>Heat treatment edu pack 1.2</h2>
        <p><b>Three integrated educational tools for steel heat treatment:</b></p>
        <ul>
            <li><b>Quenching Studio</b> – Phase transformation diagrams (TTT, CCT) and hardness prediction.</li>
            <li><b>Jomina Analyser</b> – Hardenability calculations (DI, Jominy curves) using ASTM, Just, de Cremona methods.</li>
            <li><b>Temperdata</b> – Tempering curve analysis and hardness vs. temperature plots.</li>
        </ul>
        <p><b>License:</b> GNU GPL v3</p>
        <p><b>Author:</b> Marek Góral (<a href='mailto:m_goral@interia.pl'>m_goral@interia.pl</a>)</p>
        <p><b>Free for educational and didactic use.</b></p>
        """
        QMessageBox.about(self, "About Heat treatment edu pack", about_text)

    def closeEvent(self, event):
        for w in self.windows:
            w.close()
        event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    default_font = QFont("Segoe UI", 10)
    app.setFont(default_font)
    launcher = SteelSuiteLauncher()
    launcher.show()
    sys.exit(app.exec_())
