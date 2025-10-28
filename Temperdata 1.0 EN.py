import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backend_bases import MouseEvent
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, 
                             QWidget, QComboBox, QLabel, QPushButton, QFileDialog,
                             QMessageBox, QTabWidget, QTextEdit, QTableWidget, 
                             QTableWidgetItem, QHeaderView, QSplitter)
from PyQt5.QtCore import Qt, QFile, QFileInfo
from PyQt5.QtGui import QFont
import csv
import os

"""
Temperdata - Steel Tempering Analysis Software

Copyright (C) 2025 Marek Góral

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Contact: m_goral@interia.pl
"""

class TemperingCurveApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Temperdata v 1.0 - Steel Tempering Analysis")
        self.setGeometry(100, 100, 1200, 800)
        
        # Variables for data storage
        self.data = None
        self.current_steel_type = None
        self.current_steel_data = None
        self.annotation = None
        self.plots = []
        self.plot_data = []  # Stores plot point data
        
        # Create tabs
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        # Tab 1: Plot
        self.tab1 = QWidget()
        self.tabs.addTab(self.tab1, "Plot")
        
        # Tab 2: Steel information
        self.tab2 = QWidget()
        self.tabs.addTab(self.tab2, "Steel Information")
        
        # Tab 3: Data table
        self.tab3 = QWidget()
        self.tabs.addTab(self.tab3, "Data Table")
        
        # Tab 4: Description and manual
        self.tab4 = QWidget()
        self.tabs.addTab(self.tab4, "Description & Manual")
        
        # Initialize tabs
        self.init_plot_tab()
        self.init_info_tab()
        self.init_table_tab()
        self.init_manual_tab()
        
        # Initialize empty plot
        self.initialize_plot()
        
        # Try to auto-load data
        self.try_auto_load_data()

    def init_plot_tab(self):
        """Initialize plot tab"""
        layout = QVBoxLayout(self.tab1)
        
        # Controls
        controls_layout = QHBoxLayout()
        
        self.steel_label = QLabel("Select steel grade:")
        self.steel_combo = QComboBox()
        self.steel_combo.currentTextChanged.connect(self.update_all)
        
        self.load_button = QPushButton("Load CSV Data")
        self.load_button.clicked.connect(self.load_data)
        
        self.save_plot_button = QPushButton("Save Plot")
        self.save_plot_button.clicked.connect(self.save_plot)
        self.save_plot_button.setEnabled(False)
        
        controls_layout.addWidget(self.steel_label)
        controls_layout.addWidget(self.steel_combo)
        controls_layout.addWidget(self.load_button)
        controls_layout.addWidget(self.save_plot_button)
        
        # Plot
        self.figure, self.ax = plt.subplots(figsize=(10, 6))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        
        # Add to layout
        layout.addLayout(controls_layout)
        layout.addWidget(self.canvas)

    def init_info_tab(self):
        """Initialize steel information tab"""
        layout = QVBoxLayout(self.tab2)
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setFont(QFont("Arial", 10))
        layout.addWidget(self.info_text)

    def init_table_tab(self):
        """Initialize data table tab"""
        layout = QVBoxLayout(self.tab3)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.save_table_button = QPushButton("Save Table to CSV")
        self.save_table_button.clicked.connect(self.save_table)
        self.save_table_button.setEnabled(False)
        
        button_layout.addWidget(self.save_table_button)
        button_layout.addStretch()
        
        # Table
        self.table_widget = QTableWidget()
        self.table_widget.setAlternatingRowColors(True)
        self.table_widget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        layout.addLayout(button_layout)
        layout.addWidget(self.table_widget)

    def init_manual_tab(self):
        """Initialize description and manual tab"""
        layout = QVBoxLayout(self.tab4)
        manual_text = QTextEdit()
        manual_text.setReadOnly(True)
        manual_text.setFont(QFont("Arial", 10))
        
        # Manual text
        manual_content = """
        <h1>Steel Tempering Curve Analysis - User Manual</h1>
        
        <h2>Program Description</h2>
        <p>This program is designed for visualization and analysis of steel tempering process data. 
        It enables viewing tempering curves for various steel grades, 
        reviewing their chemical composition, and exporting data and plots.</p>
        
        <h2>License Information</h2>
        <p><b>GNU GENERAL PUBLIC LICENSE Version 3</b></p>
        <p>This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.</p>
        
        <p>This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU General Public License for more details.</p>
        
        <h2>Program Features</h2>
        <ul>
            <li><b>Data Loading:</b> The program automatically attempts to load the file "Tempering data for carbon and low alloy steels - Raiipa.csv" 
            from the current directory. If the file doesn't exist, use the "Load CSV Data" button.</li>
            <li><b>Plots:</b> The "Plot" tab shows tempering curves for the selected steel grade. 
            Different line colors correspond to different tempering times.</li>
            <li><b>Tooltips:</b> Hovering the mouse over a plot point displays detailed information 
            about tempering time, temperature, and hardness.</li>
            <li><b>Steel Information:</b> The "Steel Information" tab shows the chemical composition of the selected steel grade.</li>
            <li><b>Data Table:</b> The "Data Table" tab presents all data for the selected steel in tabular form.</li>
            <li><b>Export:</b> Ability to save the plot as an image file and data as a CSV file.</li>
        </ul>
        
        <h2>User Guide</h2>
        <ol>
            <li>After starting the program, data is automatically loaded from the CSV file (if it exists in the current directory).</li>
            <li>Select a steel grade from the dropdown list.</li>
            <li>Navigate between tabs to see different data views:
                <ul>
                    <li><b>Plot:</b> Visualization of tempering curves</li>
                    <li><b>Steel Information:</b> Chemical composition of selected steel</li>
                    <li><b>Data Table:</b> Complete data in table format</li>
                </ul>
            </li>
            <li>Use the "Save Plot" and "Save Table to CSV" buttons to export data.</li>
            <li>Hover the mouse over plot points to see detailed information.</li>
        </ol>
        
        <h2>Data File Requirements</h2>
        <p>The CSV file should contain the following columns:</p>
        <ul>
            <li>Steel type - steel grade</li>
            <li>Tempering time (s) - tempering time in seconds</li>
            <li>Tempering temperature (ºC) - tempering temperature in degrees Celsius</li>
            <li>Final hardness (HRC) - post tempering - final hardness in HRC scale</li>
            <li>Chemical composition: C (%wt), Mn (%wt), P (%wt), S (%wt), Si (%wt), Ni (%wt), Cr (%wt), Mo (%wt), V (%wt), Al (%wt), Cu (%wt)</li>
        </ul>
        
        <h2>Data Source</h2>
        <p>The application is based on analysis of data collected in the database available at:</p>
        <p>https://www.kaggle.com/datasets/rgerschtzsauer/tempering-data-for-carbon-and-low-alloy-steels</p>
        
        <h2>Author</h2>
        <p>Program developed for analysis of steel tempering process data.</p>
        <p>Free software for educational use. Author is not responsible for program operation.</p>
        <p>Copyright (c) 2025 Marek Góral m_goral@interia.pl</p>
        """
        
        manual_text.setHtml(manual_content)
        layout.addWidget(manual_text)

    def initialize_plot(self):
        """Initializes empty plot"""
        self.ax.clear()
        self.ax.set_xlabel('Tempering Temperature [°C]')
        self.ax.set_ylabel('Final Hardness [HRC]')
        self.ax.set_title('Steel Tempering Curves')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        self.canvas.draw()

    def try_auto_load_data(self):
        """Tries to automatically load data from default file"""
        default_file = "Tempering data for carbon and low alloy steels - Raiipa.csv"
        
        if os.path.exists(default_file):
            try:
                self.data = pd.read_csv(default_file)
                
                # Check if required columns exist
                required_columns = ['Steel type', 'Tempering temperature (ºC)', 'Final hardness (HRC) - post tempering']
                if not all(col in self.data.columns for col in required_columns):
                    QMessageBox.critical(self, "Error", "CSV file does not contain required columns!")
                    return
                
                # Get available steel grades
                steel_types = sorted(self.data['Steel type'].unique())
                self.steel_combo.clear()
                self.steel_combo.addItems(steel_types)
                
                # Select first steel grade
                if steel_types:
                    self.steel_combo.setCurrentText(steel_types[0])
                    self.save_plot_button.setEnabled(True)
                    self.save_table_button.setEnabled(True)
                
                QMessageBox.information(self, "Success", f"Data loaded successfully from file:\n{default_file}\nNumber of rows: {len(self.data)}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load data: {str(e)}")
        else:
            QMessageBox.information(self, "Information", "Default data file not found. Please use the 'Load CSV Data' button.")

    def load_data(self):
        """Loads data from CSV file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select CSV File", "", "CSV Files (*.csv)"
        )
        
        if file_path:
            try:
                # Load data from CSV file
                self.data = pd.read_csv(file_path)
                
                # Check if required columns exist
                required_columns = ['Steel type', 'Tempering temperature (ºC)', 'Final hardness (HRC) - post tempering']
                if not all(col in self.data.columns for col in required_columns):
                    QMessageBox.critical(self, "Error", "CSV file does not contain required columns!")
                    return
                
                # Get available steel grades
                steel_types = sorted(self.data['Steel type'].unique())
                self.steel_combo.clear()
                self.steel_combo.addItems(steel_types)
                
                # Select first steel grade
                if steel_types:
                    self.steel_combo.setCurrentText(steel_types[0])
                    self.save_plot_button.setEnabled(True)
                    self.save_table_button.setEnabled(True)
                
                QMessageBox.information(self, "Success", f"Data loaded successfully!\nNumber of rows: {len(self.data)}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load data: {str(e)}")

    def update_all(self, steel_type):
        """Updates all tabs for selected steel grade"""
        if self.data is None or steel_type == "":
            return
            
        self.current_steel_type = steel_type
        
        # Filter data for selected steel grade
        self.current_steel_data = self.data[self.data['Steel type'] == steel_type]
        
        if self.current_steel_data.empty:
            return
            
        # Update plot
        self.update_plot()
        
        # Update steel information
        self.update_info()
        
        # Update table
        self.update_table()

    def update_plot(self):
        """Updates plot for selected steel grade"""
        if self.current_steel_data is None:
            return
            
        # Sort data by tempering temperature
        steel_data = self.current_steel_data.sort_values('Tempering temperature (ºC)')
        
        # Get unique tempering times
        tempering_times = steel_data['Tempering time (s)'].unique()
        
        # Prepare plot
        self.ax.clear()
        
        # Reset data storage lists
        self.plots = []
        self.plot_data = []
        
        # For each tempering time, draw a line
        for time in tempering_times:
            time_data = steel_data[steel_data['Tempering time (s)'] == time]
            plot = self.ax.plot(
                time_data['Tempering temperature (ºC)'], 
                time_data['Final hardness (HRC) - post tempering'],
                marker='o', 
                label=f'{time} s'
            )
            self.plots.append(plot[0])
            
            # Store point data for this tempering time
            for i, (temp, hardness) in enumerate(zip(time_data['Tempering temperature (ºC)'], 
                                                     time_data['Final hardness (HRC) - post tempering'])):
                self.plot_data.append({
                    'x': temp,
                    'y': hardness,
                    'time': time,
                    'plot': plot[0]
                })
        
        # Plot configuration
        self.ax.set_xlabel('Tempering Temperature [°C]')
        self.ax.set_ylabel('Final Hardness [HRC]')
        self.ax.set_title(f'Tempering Curves for {self.current_steel_type}')
        self.ax.grid(True, linestyle='--', alpha=0.7)
        self.ax.legend(title='Tempering Time')
        
        # Adjust axis scale
        self.ax.autoscale()
        
        # Refresh plot
        self.canvas.draw()

    def update_info(self):
        """Updates steel information"""
        if self.current_steel_data is None:
            return
            
        # Get first row of data for selected steel
        first_row = self.current_steel_data.iloc[0]
        
        # Prepare information text
        info_text = f"<h2>Steel Information: {self.current_steel_type}</h2>"
        info_text += "<h3>Chemical Composition:</h3>"
        info_text += "<ul>"
        
        # List of chemical elements to display
        elements = ['C (%wt)', 'Mn (%wt)', 'P (%wt)', 'S (%wt)', 'Si (%wt)', 
                   'Ni (%wt)', 'Cr (%wt)', 'Mo (%wt)', 'V (%wt)', 'Al (%wt)', 'Cu (%wt)']
        
        for element in elements:
            value = first_row.get(element, 'no data')
            if value != 'no data' and not pd.isna(value):
                element_name = element.split(' ')[0]
                info_text += f"<li>{element_name}: {value}%</li>"
        
        info_text += "</ul>"
        
        info_text += "<h3>Additional Information:</h3>"
        info_text += f"<p>Data source: {first_row.get('Source', 'no data')}</p>"
        info_text += f"<p>Number of measurements: {len(self.current_steel_data)}</p>"
        
        # Find temperature and hardness ranges
        min_temp = self.current_steel_data['Tempering temperature (ºC)'].min()
        max_temp = self.current_steel_data['Tempering temperature (ºC)'].max()
        min_hardness = self.current_steel_data['Final hardness (HRC) - post tempering'].min()
        max_hardness = self.current_steel_data['Final hardness (HRC) - post tempering'].max()
        
        info_text += f"<p>Temperature range: {min_temp} - {max_temp}°C</p>"
        info_text += f"<p>Hardness range: {min_hardness} - {max_hardness} HRC</p>"
        
        self.info_text.setHtml(info_text)

    def update_table(self):
        """Updates data table"""
        if self.current_steel_data is None:
            return
            
        # Sort data by temperature and tempering time
        sorted_data = self.current_steel_data.sort_values(['Tempering time (s)', 'Tempering temperature (ºC)'])
        
        # Prepare table
        self.table_widget.setRowCount(len(sorted_data))
        self.table_widget.setColumnCount(4)
        self.table_widget.setHorizontalHeaderLabels(['Time [s]', 'Temperature [°C]', 'Hardness [HRC]', 'Source'])
        
        # Fill table with data
        for row_idx, (_, row_data) in enumerate(sorted_data.iterrows()):
            self.table_widget.setItem(row_idx, 0, QTableWidgetItem(str(row_data['Tempering time (s)'])))
            self.table_widget.setItem(row_idx, 1, QTableWidgetItem(str(row_data['Tempering temperature (ºC)'])))
            self.table_widget.setItem(row_idx, 2, QTableWidgetItem(str(row_data['Final hardness (HRC) - post tempering'])))
            self.table_widget.setItem(row_idx, 3, QTableWidgetItem(str(row_data.get('Source', 'no data'))))

    def on_hover(self, event):
        """Handles mouse hover event on plot points"""
        if event.inaxes != self.ax or not self.plot_data:
            if self.annotation:
                self.annotation.remove()
                self.annotation = None
                self.canvas.draw()
            return
            
        # Check which point is closest to cursor
        min_distance = float('inf')
        closest_point = None
        
        for point in self.plot_data:
            distance = ((event.xdata - point['x']) ** 2 + (event.ydata - point['y']) ** 2) ** 0.5
            if distance < min_distance:
                min_distance = distance
                closest_point = point
        
        # If point found within reasonable distance (larger area)
        if min_distance < 15:  # Increased distance threshold
            # Remove existing annotation
            if self.annotation:
                self.annotation.remove()
            
            # Create new annotation
            self.annotation = self.ax.annotate(
                f'Time: {closest_point["time"]} s\nTemp: {closest_point["x"]:.1f}°C\nHardness: {closest_point["y"]:.1f} HRC',
                xy=(closest_point["x"], closest_point["y"]),
                xytext=(20, 20),
                textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', alpha=0.9),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
            )
            
            # Refresh plot
            self.canvas.draw()
        elif self.annotation:
            # If cursor is away from points, remove annotation
            self.annotation.remove()
            self.annotation = None
            self.canvas.draw()

    def save_plot(self):
        """Saves plot to file"""
        if self.current_steel_data is None:
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", f"plot_{self.current_steel_type}.png", 
            "Images (*.png *.jpg *.pdf *.svg)"
        )
        
        if file_path:
            try:
                self.figure.savefig(file_path, dpi=300, bbox_inches='tight')
                QMessageBox.information(self, "Success", f"Plot saved successfully to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save plot: {str(e)}")

    def save_table(self):
        """Saves table to CSV file"""
        if self.current_steel_data is None:
            return
            
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Table", f"table_{self.current_steel_type}.csv", 
            "CSV Files (*.csv)"
        )
        
        if file_path:
            try:
                # Sort data before saving
                sorted_data = self.current_steel_data.sort_values(['Tempering time (s)', 'Tempering temperature (ºC)'])
                
                # Save to CSV file
                sorted_data.to_csv(file_path, index=False, encoding='utf-8')
                QMessageBox.information(self, "Success", f"Table saved successfully to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save table: {str(e)}")

def main():
    app = QApplication(sys.argv)
    window = TemperingCurveApp()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
