import sys
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                             QComboBox, QLabel, QLineEdit, QPushButton, QGroupBox, QFormLayout,
                             QTabWidget, QMessageBox, QFileDialog, QTextEdit, QTableWidget, 
                             QTableWidgetItem, QHeaderView, QCheckBox, QDialog, QDialogButtonBox)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QPalette, QColor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import json
import csv
import re
from datetime import datetime

# Standard Jominy test distances (in mm) - limited to 50 mm
STANDARD_JOMINY_DISTANCES = [1.5, 3, 5, 7, 9, 11, 13, 15, 20, 25, 30, 35, 40, 45, 50]

# Hardness conversion functions
def hrc_to_hb(hrc):
    """Convert HRC to HB (approximate)"""
    if hrc < 20:
        return hrc * 10
    elif hrc > 60:
        return 650
    else:
        # Approximation formula
        return 0.363 * hrc**2 + 13.93 * hrc + 100.8

def hb_to_hrc(hb):
    """Convert HB to HRC (approximate)"""
    if hb < 200:
        return hb / 10
    elif hb > 600:
        return 65
    else:
        # Approximation formula
        return 0.094 * hb - 16.9

# Functions for calculating half-martensitic hardness according to ASTM A255
def half_martensite_hardness(carbon_content):
    """
    Calculates half-martensitic hardness based on carbon content
    According to ASTM A255, Table 7
    """
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

# ASTM/Grossman method according to algorithm from https://mxteen.github.io/metallurgical-engineers-handbook/hardenability-calculation.html
def astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu=0):
    """
    ASTM/Grossman method for calculating critical diameter DI
    Algorithm compatible with: https://mxteen.github.io/metallurgical-engineers-handbook/hardenability-calculation.html
    
    Coefficients:
    - D_I = D_I_base * f_Mn * f_Si * f_Ni * f_Cr * f_Mo (multiplicative formula)
    - D_I_base depends on carbon content
    - Alloy coefficients are multiplied by base value
    """
    
    # Determine base DI value based on carbon content
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
    
    # Calculate alloy coefficients
    f_Mn = 1.0 + 4.10 * Mn
    f_Si = 1.0 + 0.70 * Si
    f_Ni = 1.0 + 0.50 * Ni
    f_Cr = 1.0 + 2.33 * Cr
    f_Mo = 1.0 + 3.14 * Mo
    f_Cu = 1.0 + 0.52 * Cu
    
    # Calculate final DI value
    DI = D_I_base * f_Mn * f_Si * f_Ni * f_Cr * f_Mo * f_Cu
    
    # Convert from inches to mm
    return DI * 25.4

def just_DI(C, Mn, Cr, Mo, Ni, V):
    """
    Just method (Germany, 1969)
    Coefficients:
    - f_C = 0.123 + 0.289*C + 0.330*C²
    - f_Mn = 1 + 3.333*Mn
    - f_Cr = 1 + 1.167*Cr
    - f_Mo = 1 + 2.5*Mo
    - f_Ni = 1 + 0.362*Ni
    - f_V = 1 + 0.667*V
    - Constant K = 25.4 mm
    """
    K = 25.4
    f_C = 0.123 + 0.289*C + 0.330*C**2
    f_Mn = 1 + 3.333*Mn
    f_Cr = 1 + 1.167*Cr
    f_Mo = 1 + 2.5*Mo
    f_Ni = 1 + 0.362*Ni
    f_V = 1 + 0.667*V
    return K * f_C * f_Mn * f_Cr * f_Mo * f_Ni * f_V

def cremona_DI(C, Mn, Si, Cr, Ni, Mo, Cu, V, W):
    """
    de Cremona method (France, 1970)
    Coefficients:
    - f_C = 0.19 + 0.27*C
    - f_Mn = 1 + 3*Mn
    - f_Si = 1 + 0.7*Si
    - f_Cr = 1 + Cr
    - f_Ni = 1 + 0.27*Ni
    - f_Mo = 1 + 3*Mo
    - f_Cu = 1 + 0.4*Cu
    - f_V = 1 + 2*V
    - f_W = 1 + 2*W
    - Constant K = 25.4 mm
    """
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

# Function for calculating Jominy curve using Just method (simplified, without interpolation)
def just_jominy_curve(C, Mn, Si, Cr, Ni, Mo, V, grain_size, distances):
    """
    Just method for calculating Jominy curve
    According to: Journal of Materials Processing Technology 64 (1997) 117-126
    
    Coefficients:
    - C_{c1} = 88, C_{c2} = -0.0135, C_{Cr} = 19, C_{Ni} = 6.3, 
      C_{Mn} = 16, C_{Mo} = 35, C_{Si} = 5
    - m = 1/1.5875
    - N = grain size according to ASTM
    """
    # Constants
    m = 1 / 1.5875
    hardness_values = []
    
    for e in distances:
        # Main formula for all distances (without interpolation)
        term1 = 88 * C
        term2 = -0.0135 * (m * e)**2 / max(C, 0.001)  # Protection against division by zero
        term3 = 19 * Cr
        term4 = 6.3 * Ni
        term5 = 16 * Mn
        term6 = 35 * Mo
        term7 = 5 * Si
        term8 = -0.82 * grain_size
        term9 = -20 * np.sqrt(m * e)
        term10 = 2.11 * (m * e)
        constant = -2
        
        hardness = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + constant
        
        # Limit hardness to range 0-65 HRC
        hardness = max(0, min(65, hardness))
        hardness_values.append(hardness)
    
    return np.array(hardness_values)

# Function for calculating Jominy curve using ASTM A255 method
def astm_jominy_curve(C, Mn, Si, Ni, Cr, Mo, Cu, grain_size, distances):
    """
    ASTM A255 method for calculating Jominy curve
    According to: ASTM A255-02, Tables 2-18
    
    Uses calculated DI and carbon content to determine the curve
    """
    # Calculate DI using ASTM/Grossman method
    DI = astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu)
    
    # Calculate initial hardness based on carbon content (Table 7)
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
    
    # Calculate half-martensitic hardness
    half_martensitic = half_martensite_hardness(C)
    
    # Simulate Jominy curve based on DI and carbon content
    hardness_values = []
    for distance in distances:
        # Simplified model based on DI and distance
        # Higher DI means slower hardness decrease with distance
        # Higher carbon content means higher initial hardness
        decay_factor = np.exp(-distance / (DI / 10))
        hardness = half_martensitic + (initial_hardness - half_martensitic) * decay_factor
        hardness = max(20, min(65, hardness))  # Limit to range 20-65 HRC
        hardness_values.append(hardness)
    
    return np.array(hardness_values)

# Function for calculating ideal critical diameter from actual
def calculate_ideal_DI(actual_DI, severity_factor=1.0):
    """
    Calculates ideal critical diameter based on actual critical diameter
    and cooling intensity factor (severity factor)
    
    D_ideal = D_actual * (H / 0.8)^(1/2)
    where H is cooling intensity factor (for ideal cooling H=∞, but we assume H=1.0 for water)
    """
    return actual_DI * (1.0 / 0.8) ** 0.5

# Function for converting Jominy distance to DI according to ASTM A255
def jominy_distance_to_DI(j_distance):
    """
    Converts Jominy distance to critical diameter DI according to ASTM A255
    Uses data from Table 9 (mm) of ASTM A255 standard
    """
    # Data from Table 9 ASTM A255 (mm)
    jominy_data = {
        'distance': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
                     16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
                     31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0,
                     46.0, 47.0, 48.0, 49.0, 50.0],
        'DI': [8.4, 15.7, 22.9, 29.7, 36.3, 42.9, 48.2, 54.2, 59.5, 64.2, 68.6, 72.1, 76.4, 80.1, 84.0,
               87.6, 90.1, 94.2, 97.1, 100.5, 103.7, 106.5, 109.7, 112.2, 114.9, 117.4, 119.9, 122.4, 124.7, 127.1,
               129.0, 131.4, 133.5, 135.2, 137.1, 139.1, 140.9, 142.8, 144.7, 146.4, 148.3, 150.1, 151.7, 153.4, 154.1,
               156.5, 157.8, 159.2, 160.5, 161.8]
    }
    
    # Linear interpolation
    return np.interp(j_distance, jominy_data['distance'], jominy_data['DI'])

# Steel database with descriptions, heat treatment conditions and hardness after quenching
steel_database = {
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
        "description": "Case hardening and heat treatable steel, used for machine components, gears",
        "heat_treatment": "Quenching: 830-850°C in oil, Tempering: 550-650°C",
        "hardness": "200-250 HB (case hardened), 250-300 HB (heat treated)",
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
        "description": "Case hardening steel, used for gears, shafts, transmission components",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (carburized surface), 280-320 HB (core)",
        "hardness_after_quenching": "58-62 HRC (surface)"
    },
    "1.7131; 16MnCr5; 16HG": {
        "C": 0.16, "Mn": 1.10, "Si": 0.30, "Cr": 1.00, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 70, "half_martensitic": 30,
        "description": "Case hardening steel, used for gears, shafts, machine components",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (carburized surface), 250-300 HB (core)",
        "hardness_after_quenching": "58-62 HRC (surface)"
    },
    "1.7147; 20MnCr5; 20HG": {
        "C": 0.20, "Mn": 1.10, "Si": 0.30, "Cr": 1.10, "Ni": 0.10, "Mo": 0.05, "Cu": 0.10, "V": 0.0, "W": 0.0,
        "grain_size": 7, "typical_DI": 75, "half_martensitic": 32,
        "description": "Case hardening steel, used for gears, shafts, machine components",
        "heat_treatment": "Carburizing: 900-930°C, Quenching: 820-850°C in oil, Tempering: 180-200°C",
        "hardness": "58-62 HRC (carburized surface), 280-320 HB (core)",
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

class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("About Jomina Analyzer")
        self.setModal(True)
        self.resize(500, 400)
        
        layout = QVBoxLayout()
        
        # Title
        title = QLabel("Jomina Analyzer")
        title.setStyleSheet("font-size: 18pt; font-weight: bold; color: #2E86AB;")
        layout.addWidget(title)
        
        # Version
        version = QLabel("Version 1.0")
        version.setStyleSheet("font-size: 12pt; color: #666;")
        layout.addWidget(version)
        
        # Description
        description = QTextEdit()
        description.setReadOnly(True)
        description.setPlainText("""
Jominy Analyzer is a comprehensive software tool for metallurgical engineers and heat treatment specialists.

MAIN FEATURES:
• Calculation of critical diameter (DI) using multiple methods (ASTM/Grossman, Just, de Cremona)
• Prediction of Jominy hardenability curves based on chemical composition
• Comparison of calculated curves with experimental data
• Support for up to 4 experimental hardness curves
• Steel database with comprehensive material information
• Professional visualization of results

CALCULATION METHODS:
• ASTM/Grossman method for DI calculation
• Just method (Germany, 1969) for DI and Jominy curve
• de Cremona method (France, 1970) for DI
• ASTM A255 standard for Jominy curve prediction

This software is developed under GNU General Public License v3.0
        """)
        layout.addWidget(description)
        
        # License info
        license_info = QLabel("License: GNU General Public License v3.0")
        license_info.setStyleSheet("font-weight: bold; color: #A23B72;")
        layout.addWidget(license_info)
        
        # Buttons
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
        self.axes.set_xlim(0, 50)  # Limit X-axis to 50 mm
        self.curves = []  # Stores curve information
        
    def plot_jominy(self, distance_just, hardness_just, distance_astm, hardness_astm, 
                   experimental_data=None, steel_name="", half_martensitic=None, show_just=True):
        # Clear all curves
        self.curves = []
        
        # Add Just curve if enabled
        if show_just:
            self.curves.append({
                'distance': distance_just,
                'hardness': hardness_just,
                'steel_name': steel_name,
                'curve_name': "Just Method",
                'color': 'r',
                'linestyle': '-'
            })
        
        # Add ASTM curve
        self.curves.append({
            'distance': distance_astm,
            'hardness': hardness_astm,
            'steel_name': steel_name,
            'curve_name': "ASTM Method",
            'color': 'b',
            'linestyle': '--'
        })
        
        # Add experimental data if available
        if experimental_data:
            colors = ['g', 'm', 'c', 'y']
            markers = ['o', 's', '^', 'D']
            for i, (distance_exp, hardness_exp) in enumerate(experimental_data):
                if distance_exp is not None and hardness_exp is not None and len(distance_exp) > 0:
                    self.curves.append({
                        'distance': distance_exp,
                        'hardness': hardness_exp,
                        'steel_name': steel_name,
                        'curve_name': f"Experimental Curve {i+1}",
                        'color': colors[i % len(colors)],
                        'linestyle': '-',
                        'marker': markers[i % len(markers)]
                    })
            
        self.update_plot(half_martensitic)
    
    def update_plot(self, half_martensitic=None):
        self.axes.clear()
        
        for curve in self.curves:
            if 'marker' in curve:
                # For experimental data - line with markers
                self.axes.plot(curve['distance'], curve['hardness'], 
                              color=curve['color'], linestyle=curve['linestyle'], 
                              marker=curve['marker'], markersize=4, linewidth=2,
                              label=f"{curve['curve_name']}")
            else:
                # For calculated curves - line only
                self.axes.plot(curve['distance'], curve['hardness'], 
                              color=curve['color'], linestyle=curve['linestyle'], linewidth=2,
                              label=f"{curve['curve_name']}")
        
        # Add horizontal line for half-martensitic hardness
        if half_martensitic is not None:
            self.axes.axhline(y=half_martensitic, color='g', linestyle='--', linewidth=2, 
                             label=f'Half-martensitic hardness: {half_martensitic} HRC')
        
        self.axes.grid(True)
        self.axes.set_xlabel('Distance from quenched end [mm]')
        self.axes.set_ylabel('Hardness [HRC]')
        self.axes.set_title('Jominy Curves Comparison')
        self.axes.set_xlim(0, 50)  # Limit X-axis to 50 mm
        self.axes.legend()
        self.draw()
    
    def clear_curves(self):
        self.curves = []
        self.update_plot()

class DIComparisonPlot(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        self.axes = self.fig.add_subplot(111)
        
    def plot_comparison(self, methods, di_values, steel_name, jominy_di=None, ideal_di=None):
        self.axes.clear()
        y_pos = np.arange(len(methods))
        bars = self.axes.barh(y_pos, di_values, align='center', alpha=0.7)
        self.axes.set_yticks(y_pos)
        self.axes.set_yticklabels(methods)
        self.axes.set_xlabel('Critical Diameter DI [mm]')
        self.axes.set_title(f'DI Calculation Methods Comparison for {steel_name}')
        
        # Add values on bars
        for i, v in enumerate(di_values):
            self.axes.text(v + 1, i, f'{v:.2f}', va='center')
        
        # Add line for Jominy curve value if available
        if jominy_di is not None:
            self.axes.axvline(x=jominy_di, color='r', linestyle='--', linewidth=2, label=f'Actual DI: {jominy_di:.2f} mm')
        
        # Add line for ideal DI value if available
        if ideal_di is not None:
            self.axes.axvline(x=ideal_di, color='g', linestyle='--', linewidth=2, label=f'Ideal DI: {ideal_di:.2f} mm')
        
        if jominy_di is not None or ideal_di is not None:
            self.axes.legend()
        
        self.draw()

class SteelCalculatorApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Steel Hardenability Calculator - Jominy Analyzer")
        self.setGeometry(100, 100, 1600, 900)  # Larger window for Full HD
        
        # Apply professional style
        self.apply_style()
        
        # Data for storing Jominy curves
        self.jominy_curves = []  # List to store multiple curves
        self.jominy_di = None
        self.ideal_di = None
        self.show_just_method = True  # Default show Just method
        
        # Store all calculated DI values
        self.di_results = {}
        
        # Central widget with tabs
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QVBoxLayout(central_widget)
        
        # Create tabs
        self.tabs = QTabWidget()
        self.tab1 = QWidget()  # DI Calculations
        self.tab2 = QWidget()  # Jominy Curve
        self.tab3 = QWidget()  # DI Methods Comparison
        self.tab4 = QWidget()  # Methods Description
        self.tab5 = QWidget()  # About
        
        self.tabs.addTab(self.tab1, "DI Calculations")
        self.tabs.addTab(self.tab2, "Jominy Curve")
        self.tabs.addTab(self.tab3, "DI Methods Comparison")
        self.tabs.addTab(self.tab4, "Methods Description")
        self.tabs.addTab(self.tab5, "About")
        
        # Setup each tab
        self.setup_tab1()
        self.setup_tab2()
        self.setup_tab3()
        self.setup_tab4()
        self.setup_tab5()
        
        main_layout.addWidget(self.tabs)
        
        # Initialize with first steel
        self.on_steel_changed(self.steel_combo.currentText())
    
    def apply_style(self):
        """Apply professional styling to the application"""
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f0f0f0;
            }
            QWidget {
                font-family: Segoe UI, Arial, sans-serif;
                font-size: 9pt;
            }
            QGroupBox {
                font-weight: bold;
                border: 1px solid #cccccc;
                border-radius: 5px;
                margin-top: 1ex;
                padding-top: 10px;
                background-color: white;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px 0 5px;
                color: #2E86AB;
            }
            QPushButton {
                background-color: #4CAF50;
                border: none;
                color: white;
                padding: 6px 12px;
                text-align: center;
                text-decoration: none;
                font-size: 9pt;
                margin: 2px 1px;
                border-radius: 4px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QPushButton:pressed {
                background-color: #3d8b40;
            }
            QLineEdit, QComboBox {
                padding: 4px;
                border: 1px solid #cccccc;
                border-radius: 3px;
                background-color: white;
                font-size: 9pt;
            }
            QTabWidget::pane {
                border: 1px solid #cccccc;
                background-color: white;
            }
            QTabBar::tab {
                background-color: #e0e0e0;
                padding: 8px 16px;
                margin-right: 2px;
                border-top-left-radius: 4px;
                border-top-right-radius: 4px;
                font-size: 9pt;
            }
            QTabBar::tab:selected {
                background-color: white;
                font-weight: bold;
                color: #2E86AB;
            }
            QTableWidget {
                gridline-color: #d0d0d0;
                selection-background-color: #4CAF50;
                font-size: 9pt;
            }
            QHeaderView::section {
                background-color: #e0e0e0;
                padding: 4px;
                border: 1px solid #cccccc;
                font-weight: bold;
            }
            QTextEdit {
                border: 1px solid #cccccc;
                border-radius: 3px;
                background-color: white;
                font-size: 9pt;
            }
            QLabel {
                font-size: 9pt;
            }
        """)
    
    def setup_tab1(self):
        # Tab 1 layout - Three columns for better space utilization
        layout = QHBoxLayout(self.tab1)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # Left column - Input controls
        left_column = QVBoxLayout()
        left_column.setSpacing(10)
        
        # Steel selection
        steel_group = QGroupBox("Steel Selection")
        steel_layout = QVBoxLayout(steel_group)
        
        self.steel_combo = QComboBox()
        self.steel_combo.addItems(steel_database.keys())
        self.steel_combo.currentTextChanged.connect(self.on_steel_changed)
        steel_layout.addWidget(QLabel("Select steel:"))
        steel_layout.addWidget(self.steel_combo)
        
        # Chemical composition inputs
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
        
        comp_layout.addRow("Carbon (C):", self.c_input)
        comp_layout.addRow("Manganese (Mn):", self.mn_input)
        comp_layout.addRow("Silicon (Si):", self.si_input)
        comp_layout.addRow("Chromium (Cr):", self.cr_input)
        comp_layout.addRow("Nickel (Ni):", self.ni_input)
        comp_layout.addRow("Molybdenum (Mo):", self.mo_input)
        comp_layout.addRow("Copper (Cu):", self.cu_input)
        comp_layout.addRow("Vanadium (V):", self.v_input)
        comp_layout.addRow("Tungsten (W):", self.w_input)
        
        # Grain size and hardness after quenching inputs
        params_group = QGroupBox("Parameters")
        params_layout = QFormLayout(params_group)
        
        self.grain_size_input = QLineEdit("7")
        self.hardness_after_quenching_input = QLineEdit("50")
        
        params_layout.addRow("ASTM Grain Size:", self.grain_size_input)
        params_layout.addRow("Hardness after quenching [HRC]:", self.hardness_after_quenching_input)
        
        # Buttons
        button_layout = QHBoxLayout()
        self.calc_button = QPushButton("Calculate DI")
        self.calc_button.clicked.connect(self.calculate_DI)
        self.save_comp_button = QPushButton("Save Composition")
        self.save_comp_button.clicked.connect(self.save_composition)
        self.load_comp_button = QPushButton("Load Composition")
        self.load_comp_button.clicked.connect(self.load_composition)
        
        button_layout.addWidget(self.calc_button)
        button_layout.addWidget(self.save_comp_button)
        button_layout.addWidget(self.load_comp_button)
        
        # Add widgets to left column
        left_column.addWidget(steel_group)
        left_column.addWidget(comp_group)
        left_column.addWidget(params_group)
        left_column.addLayout(button_layout)
        left_column.addStretch()
        
        # Middle column - Steel information
        middle_column = QVBoxLayout()
        middle_column.setSpacing(10)
        
        info_group = QGroupBox("Steel Information")
        info_layout = QVBoxLayout(info_group)
        
        self.steel_description = QTextEdit()
        self.steel_description.setReadOnly(True)
        self.steel_description.setMaximumHeight(200)
        info_layout.addWidget(self.steel_description)
        
        # Heat treatment info
        heat_treatment_group = QGroupBox("Heat Treatment and Hardness")
        heat_treatment_layout = QVBoxLayout(heat_treatment_group)
        
        self.heat_treatment_info = QTextEdit()
        self.heat_treatment_info.setReadOnly(True)
        self.heat_treatment_info.setMaximumHeight(200)
        heat_treatment_layout.addWidget(self.heat_treatment_info)
        
        # Add to middle column
        middle_column.addWidget(info_group)
        middle_column.addWidget(heat_treatment_group)
        middle_column.addStretch()
        
        # Right column - Results
        right_column = QVBoxLayout()
        right_column.setSpacing(10)
        
        # DI results
        results_group = QGroupBox("DI Calculation Results [mm]")
        results_layout = QFormLayout(results_group)
        
        self.astm_grossman_result = QLabel("—")
        self.just_result = QLabel("—")
        self.cremona_result = QLabel("—")
        self.experimental_result = QLabel("—")
        self.ideal_experimental_result = QLabel("—")
        
        results_layout.addRow("ASTM/Grossman:", self.astm_grossman_result)
        results_layout.addRow("Just (1969):", self.just_result)
        results_layout.addRow("de Cremona (1970):", self.cremona_result)
        results_layout.addRow("Experimental (actual):", self.experimental_result)
        results_layout.addRow("Experimental (ideal):", self.ideal_experimental_result)
        
        # Save results button
        self.save_results_button = QPushButton("Save Results")
        self.save_results_button.clicked.connect(self.save_results)
        
        # Add to right column
        right_column.addWidget(results_group)
        right_column.addWidget(self.save_results_button)
        right_column.addStretch()
        
        # Set column proportions
        left_column_widget = QWidget()
        left_column_widget.setLayout(left_column)
        left_column_widget.setMaximumWidth(400)
        
        middle_column_widget = QWidget()
        middle_column_widget.setLayout(middle_column)
        middle_column_widget.setMaximumWidth(350)
        
        right_column_widget = QWidget()
        right_column_widget.setLayout(right_column)
        right_column_widget.setMaximumWidth(300)
        
        # Add to main layout
        layout.addWidget(left_column_widget)
        layout.addWidget(middle_column_widget)
        layout.addWidget(right_column_widget)
    
    def setup_tab2(self):
        # Tab 2 layout
        layout = QHBoxLayout(self.tab2)
        layout.setSpacing(10)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # Left panel with controls - 50% width
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_panel.setMaximumWidth(800)  # 50% of 1600px
        
        # Jominy data inputs
        jominy_group = QGroupBox("Jominy Curve Data")
        jominy_layout = QVBoxLayout(jominy_group)
        
        # Curve selection checkboxes
        curve_selection_layout = QHBoxLayout()
        self.curve_checkboxes = []
        for i in range(4):
            checkbox = QCheckBox(f"Curve {i+1}")
            checkbox.setChecked(i == 0)  # First curve enabled by default
            checkbox.toggled.connect(self.update_jominy_plot)
            self.curve_checkboxes.append(checkbox)
            curve_selection_layout.addWidget(checkbox)
        
        # Show Just method checkbox
        self.show_just_checkbox = QCheckBox("Show Just Method")
        self.show_just_checkbox.setChecked(True)
        self.show_just_checkbox.toggled.connect(self.toggle_just_method)
        curve_selection_layout.addWidget(self.show_just_checkbox)
        
        jominy_layout.addLayout(curve_selection_layout)
        
        # Table for Jominy data
        self.jominy_table = QTableWidget()
        self.jominy_table.setColumnCount(5)  # Distance + 4 hardness curves
        self.jominy_table.setHorizontalHeaderLabels(["Distance [mm]", "Hardness 1 [HRC]", "Hardness 2 [HRC]", "Hardness 3 [HRC]", "Hardness 4 [HRC]"])
        self.jominy_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.jominy_table.setMinimumHeight(400)
        
        # Buttons for Jominy data
        button_layout = QHBoxLayout()
        self.add_row_button = QPushButton("Add Row")
        self.add_row_button.clicked.connect(self.add_jominy_row)
        self.remove_row_button = QPushButton("Remove Row")
        self.remove_row_button.clicked.connect(self.remove_jominy_row)
        self.calc_di_button = QPushButton("Calculate DI from Curve")
        self.calc_di_button.clicked.connect(self.calculate_DI_from_jominy)
        
        button_layout.addWidget(self.add_row_button)
        button_layout.addWidget(self.remove_row_button)
        button_layout.addWidget(self.calc_di_button)
        
        jominy_layout.addWidget(self.jominy_table)
        jominy_layout.addLayout(button_layout)
        
        # File operations
        file_group = QGroupBox("File Operations")
        file_layout = QHBoxLayout(file_group)
        
        self.load_jominy_button = QPushButton("Load Jominy Data")
        self.load_jominy_button.clicked.connect(self.load_jominy_data)
        self.save_jominy_button = QPushButton("Save Jominy Data")
        self.save_jominy_button.clicked.connect(self.save_jominy_data)
        self.save_plot_button = QPushButton("Save Plot")
        self.save_plot_button.clicked.connect(self.save_jominy_plot)
        
        file_layout.addWidget(self.load_jominy_button)
        file_layout.addWidget(self.save_jominy_button)
        file_layout.addWidget(self.save_plot_button)
        
        # Add widgets to left layout
        left_layout.addWidget(jominy_group)
        left_layout.addWidget(file_group)
        left_layout.addStretch()
        
        # Right panel with plot - 50% width
        self.jominy_plot = JominyPlot(self, width=8, height=8, dpi=100)
        
        # Add to main layout
        layout.addWidget(left_panel)
        layout.addWidget(self.jominy_plot)
        
        # Initialize table with standard distances
        self.fill_standard_distances()
        
        # Connect table changes to update plot
        self.jominy_table.cellChanged.connect(self.update_jominy_plot)
    
    def setup_tab3(self):
        # Tab 3 layout
        layout = QVBoxLayout(self.tab3)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # DI comparison plot
        self.di_comparison_plot = DIComparisonPlot(self, width=12, height=8, dpi=100)
        layout.addWidget(self.di_comparison_plot)
        
        # Button to update comparison
        self.update_comparison_button = QPushButton("Update Comparison")
        self.update_comparison_button.clicked.connect(self.update_di_comparison)
        layout.addWidget(self.update_comparison_button)
    
    def setup_tab4(self):
        # Tab 4 layout
        layout = QVBoxLayout(self.tab4)
        layout.setContentsMargins(10, 10, 10, 10)
        
        # Method description text
        method_text = QTextEdit()
        method_text.setReadOnly(True)
        
        # Load method description from file or set directly
        method_description = """
DESCRIPTION OF CRITICAL HARDENING DIAMETER (DI) CALCULATION METHODS

1. ASTM/Grossman Method
Algorithm compatible with: https://mxteen.github.io/metallurgical-engineers-handbook/hardenability-calculation.html
Formula: DI = D_I_base * f_Mn * f_Si * f_Ni * f_Cr * f_Mo * f_Cu
Coefficients:
- D_I_base depends on carbon content (0.05-0.38 inch)
- f_Mn = 1.0 + 4.10 * Mn
- f_Si = 1.0 + 0.70 * Si
- f_Ni = 1.0 + 0.50 * Ni
- f_Cr = 1.0 + 2.33 * Cr
- f_Mo = 1.0 + 3.14 * Mo
- f_Cu = 1.0 + 0.52 * Cu
- Result converted from inches to mm (× 25.4)

2. Just Method (Germany, 1969)
Formula: DI = K * f_C * f_Mn * f_Cr * f_Mo * f_Ni * f_V
Coefficients:
- K = 25.4 mm
- f_C = 0.123 + 0.289*C + 0.330*C²
- f_Mn = 1 + 3.333*Mn
- f_Cr = 1 + 1.167*Cr
- f_Mo = 1 + 2.5*Mo
- f_Ni = 1 + 0.362*Ni
- f_V = 1 + 0.667*V

3. de Cremona Method (France, 1970)
Formula: DI = K * f_C * f_Mn * f_Si * f_Cr * f_Ni * f_Mo * f_Cu * f_V * f_W
Coefficients:
- K = 25.4 mm
- f_C = 0.19 + 0.27*C
- f_Mn = 1 + 3*Mn
- f_Si = 1 + 0.7*Si
- f_Cr = 1 + Cr
- f_Ni = 1 + 0.27*Ni
- f_Mo = 1 + 3*Mo
- f_Cu = 1 + 0.4*Cu
- f_V = 1 + 2*V
- f_W = 1 + 2*W

NOTE: All chemical components are given in weight %.

JUST METHOD FOR CALCULATING JOMINY CURVE:
Just method (1969) allows direct calculation of Jominy curve based on steel chemical composition.
Formula: H_e = C_{c1}*C + C_{c2}*(m*e)^2/C + C_{Cr}*Cr + C_{Ni}*Ni + C_{Mn}*Mn + C_{Mo}*Mo + C_{Si}*Si - 0.82*N - 20*sqrt(m*e) + 2.11*(m*e) - 2
Where:
- H_e: HRC hardness at distance e from sample end
- e: distance from sample end [mm]
- m: constant = 1/1.5875
- N: ASTM grain size
- Coefficients: C_{c1}=88, C_{c2}=-0.0135, C_{Cr}=19, C_{Ni}=6.3, C_{Mn}=16, C_{Mo}=35, C_{Si}=5

ASTM A255 METHOD FOR CALCULATING JOMINY CURVE:
ASTM A255 method uses calculated DI and carbon content to determine Jominy curve.
Initial hardness is determined based on carbon content, then hardness at different distances
is calculated based on dividing coefficients dependent on DI.

CALCULATING IDEAL CRITICAL DIAMETER:
Ideal critical diameter (D_ideal) is calculated based on actual critical diameter (D_actual)
and cooling intensity factor (H) according to formula:
D_ideal = D_actual * (H / 0.8)^(1/2)
For ideal cooling, H=1.0 is assumed.
        """
        
        method_text.setPlainText(method_description)
        layout.addWidget(method_text)
    
    def setup_tab5(self):
        # Tab 5 layout - About
        layout = QVBoxLayout(self.tab5)
        layout.setContentsMargins(20, 20, 20, 20)
        
        # About information
        about_text = QTextEdit()
        about_text.setReadOnly(True)
        about_text.setHtml("""
        <html>
        <head>
        <style>
        h1 { color: #2E86AB; font-size: 24px; }
        h2 { color: #A23B72; font-size: 18px; }
        p { font-size: 11pt; line-height: 1.5; }
        li { margin-bottom: 8px; }
        .license { background-color: #f8f8f8; padding: 10px; border-radius: 5px; }
        </style>
        </head>
        <body>
        <h1>Jominy Analyzer</h1>
        <h2>Professional Steel Hardenability Analysis Tool</h2>
        
        <p><strong>Version:</strong> 1.0</p>
        <p><strong>Developed for:</strong> Metallurgical engineers, heat treatment specialists, and materials scientists</p>
        
        <h2>Main Features</h2>
        <ul>
        <li><strong>Critical Diameter (DI) Calculation:</strong> Multiple calculation methods including ASTM/Grossman, Just (1969), and de Cremona (1970)</li>
        <li><strong>Jominy Curve Prediction:</strong> Calculate hardenability curves based on chemical composition</li>
        <li><strong>Experimental Data Comparison:</strong> Support for up to 4 experimental hardness curves with visualization</li>
        <li><strong>Steel Database:</strong> Comprehensive database of common steel grades with composition and heat treatment data</li>
        <li><strong>Professional Visualization:</strong> High-quality plots for Jominy curves and DI comparison</li>
        <li><strong>Data Export:</strong> Save results in multiple formats (CSV, JSON, images)</li>
        </ul>
        
        <h2>Calculation Methods</h2>
        <ul>
        <li><strong>ASTM/Grossman:</strong> Industry standard for DI calculation based on carbon content and alloying elements</li>
        <li><strong>Just Method:</strong> German method (1969) for both DI and Jominy curve calculation</li>
        <li><strong>de Cremona Method:</strong> French method (1970) for DI calculation with comprehensive alloy factors</li>
        <li><strong>ASTM A255:</strong> Standard method for Jominy curve prediction</li>
        </ul>
        
        <div class="license">
        <h2>License Information</h2>
        <p><strong>GNU General Public License v3.0</strong></p>
        <p>This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.</p>
        
        <p>This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.</p>
        
        <p>You should have received a copy of the GNU General Public License along with this program. If not, see <a href="https://www.gnu.org/licenses/">https://www.gnu.org/licenses/</a>.</p>
        
        <p><strong>Permissions:</strong></p>
        <ul>
        <li>Commercial use</li>
        <li>Modification</li>
        <li>Distribution</li>
        <li>Patent use</li>
        <li>Private use</li>
        </ul>
        
        <p><strong>Conditions:</strong></p>
        <ul>
        <li>License and copyright notice</li>
        <li>State changes</li>
        <li>Disclose source</li>
        <li>Same license</li>
        </ul>
        </div>
        
        <h2>Technical Specifications</h2>
        <ul>
        <li><strong>Programming Language:</strong> Python 3.x</li>
        <li><strong>GUI Framework:</strong> PyQt5</li>
        <li><strong>Plotting Library:</strong> Matplotlib</li>
        <li><strong>Numerical Computing:</strong> NumPy, Pandas</li>
        <li><strong>Recommended Resolution:</strong> Full HD (1920x1080) or higher</li>
        </ul>
        
        <p><em>For technical support or feature requests, please contact the development team.</em></p>
        </body>
        </html>
        """)
        
        layout.addWidget(about_text)
    
    def toggle_just_method(self, checked):
        """Toggle display of Just method on Jominy plot"""
        self.show_just_method = checked
        self.update_jominy_plot()
    
    def fill_standard_distances(self):
        """Fills table with standard Jominy test distances"""
        self.jominy_table.setRowCount(len(STANDARD_JOMINY_DISTANCES))
        for i, distance in enumerate(STANDARD_JOMINY_DISTANCES):
            self.jominy_table.setItem(i, 0, QTableWidgetItem(str(distance)))
            for j in range(1, 5):
                self.jominy_table.setItem(i, j, QTableWidgetItem(""))
    
    def extract_hardness_value(self, hardness_str):
        """Extracts hardness value from text (e.g., '50-55 HRC' -> 52.5)"""
        if not hardness_str or "Not applicable" in hardness_str:
            return None
        
        # Find all numbers in text
        numbers = re.findall(r'\d+\.?\d*', hardness_str)
        if numbers:
            # Convert to float and calculate average
            numbers = [float(num) for num in numbers]
            return sum(numbers) / len(numbers)
        return None
    
    def on_steel_changed(self, steel_name):
        if steel_name in steel_database:
            composition = steel_database[steel_name]
            self.c_input.setText(str(composition["C"]))
            self.mn_input.setText(str(composition["Mn"]))
            self.si_input.setText(str(composition["Si"]))
            self.cr_input.setText(str(composition["Cr"]))
            self.ni_input.setText(str(composition["Ni"]))
            self.mo_input.setText(str(composition["Mo"]))
            self.cu_input.setText(str(composition["Cu"]))
            self.v_input.setText(str(composition["V"]))
            self.w_input.setText(str(composition["W"]))
            self.grain_size_input.setText(str(composition.get("grain_size", 7)))
            
            # Set hardness after quenching from database
            hardness_after_quenching = self.extract_hardness_value(composition.get("hardness_after_quenching", "50"))
            if hardness_after_quenching is not None:
                self.hardness_after_quenching_input.setText(str(hardness_after_quenching))
            
            # Update steel description
            description = f"Description: {composition['description']}\n\nTypical DI value: {composition.get('typical_DI', 'N/A')} mm"
            self.steel_description.setPlainText(description)
            
            # Update heat treatment information
            heat_treatment_info = f"Heat treatment conditions:\n{composition['heat_treatment']}\n\nHardness after heat treatment:\n{composition['hardness']}\n\nHardness after quenching:\n{composition['hardness_after_quenching']}"
            self.heat_treatment_info.setPlainText(heat_treatment_info)
            
            # Update Jominy plot
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
            grain_size = float(self.grain_size_input.text())
            hardness_after_quenching = float(self.hardness_after_quenching_input.text())
            return C, Mn, Si, Cr, Ni, Mo, Cu, V, W, grain_size, hardness_after_quenching
        except ValueError:
            QMessageBox.warning(self, "Error", "Invalid input data. Check chemical composition values.")
            return None
    
    def calculate_DI(self):
        composition = self.get_composition()
        if composition is None:
            return
        
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, grain_size, hardness_after_quenching = composition
        
        # Calculate DI using different methods
        di_astm_grossman = astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu)
        di_just = just_DI(C, Mn, Cr, Mo, Ni, V)
        di_cremona = cremona_DI(C, Mn, Si, Cr, Ni, Mo, Cu, V, W)
        
        # Save results
        self.di_results = {
            "ASTM/Grossman": di_astm_grossman,
            "Just": di_just,
            "de Cremona": di_cremona,
            "Steel": self.steel_combo.currentText()
        }
        
        # Display results
        self.astm_grossman_result.setText(f"{di_astm_grossman:.2f}")
        self.just_result.setText(f"{di_just:.2f}")
        self.cremona_result.setText(f"{di_cremona:.2f}")
        
        # Update comparison plot
        self.update_di_comparison()
        
        # Update Jominy plot
        self.update_jominy_plot()
    
    def update_di_comparison(self):
        composition = self.get_composition()
        if composition is None:
            return
        
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, grain_size, hardness_after_quenching = composition
        
        # Calculate DI using different methods
        di_astm_grossman = astm_grossman_DI(C, Mn, Si, Ni, Cr, Mo, Cu)
        di_just = just_DI(C, Mn, Cr, Mo, Ni, V)
        di_cremona = cremona_DI(C, Mn, Si, Cr, Ni, Mo, Cu, V, W)
        
        methods = ["ASTM/Grossman", "Just", "de Cremona"]
        di_values = [di_astm_grossman, di_just, di_cremona]
        
        steel_name = self.steel_combo.currentText()
        self.di_comparison_plot.plot_comparison(methods, di_values, steel_name, self.jominy_di, self.ideal_di)
    
    def add_jominy_row(self):
        row_count = self.jominy_table.rowCount()
        self.jominy_table.insertRow(row_count)
        for i in range(5):
            self.jominy_table.setItem(row_count, i, QTableWidgetItem(""))
    
    def remove_jominy_row(self):
        current_row = self.jominy_table.currentRow()
        if current_row >= 0:
            self.jominy_table.removeRow(current_row)
    
    def get_jominy_data(self):
        """Returns list of experimental curves (distance, hardness) for enabled curves"""
        distance = []
        hardness_curves = [[] for _ in range(4)]
        
        for row in range(self.jominy_table.rowCount()):
            dist_item = self.jominy_table.item(row, 0)
            
            if dist_item and dist_item.text():
                try:
                    dist_val = float(dist_item.text())
                    distance.append(dist_val)
                    
                    for i in range(4):
                        hard_item = self.jominy_table.item(row, i+1)
                        if hard_item and hard_item.text():
                            hard_val = float(hard_item.text())
                            hardness_curves[i].append(hard_val)
                        else:
                            hardness_curves[i].append(np.nan)
                except ValueError:
                    continue
        
        # Filter enabled curves
        enabled_curves = []
        for i in range(4):
            if self.curve_checkboxes[i].isChecked():
                # Remove rows with NaN values
                valid_indices = [j for j, val in enumerate(hardness_curves[i]) if not np.isnan(val)]
                if valid_indices:
                    curve_distance = np.array([distance[j] for j in valid_indices])
                    curve_hardness = np.array([hardness_curves[i][j] for j in valid_indices])
                    enabled_curves.append((curve_distance, curve_hardness))
        
        return enabled_curves
    
    def calculate_DI_from_jominy(self):
        experimental_curves = self.get_jominy_data()
        
        if not experimental_curves:
            QMessageBox.warning(self, "Error", "No experimental data available for DI calculation.")
            return
        
        composition = self.get_composition()
        if composition is None:
            return
        
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, grain_size, hardness_after_quenching = composition
        
        # Use first enabled curve for DI calculation
        distance, hardness = experimental_curves[0]
        
        if len(distance) < 2:
            QMessageBox.warning(self, "Error", "Insufficient data to calculate DI.")
            return
        
        # Sort data by distance
        sorted_indices = np.argsort(distance)
        distance = distance[sorted_indices]
        hardness = hardness[sorted_indices]
        
        # Find distance where hardness drops to half-martensitic hardness
        target_hardness = hardness_after_quenching
        
        # Interpolate to find actual DI
        if hardness[0] < target_hardness:
            # If first value is already below target hardness
            di_value = 0
        else:
            # Find point where hardness drops below target value
            idx = np.where(hardness < target_hardness)[0]
            if len(idx) > 0:
                first_below = idx[0]
                if first_below == 0:
                    di_value = 0
                else:
                    # Linear interpolation
                    x1 = distance[first_below - 1]
                    x2 = distance[first_below]
                    y1 = hardness[first_below - 1]
                    y2 = hardness[first_below]
                    
                    di_value = x1 + (x2 - x1) * (target_hardness - y1) / (y2 - y1)
            else:
                # All values are above target hardness
                di_value = distance[-1]
        
        # Convert Jominy distance to DI according to ASTM A255
        di_value = jominy_distance_to_DI(di_value)
        
        # Calculate ideal DI
        ideal_di_value = calculate_ideal_DI(di_value)
        
        self.jominy_di = di_value
        self.ideal_di = ideal_di_value
        self.experimental_result.setText(f"{di_value:.2f}")
        self.ideal_experimental_result.setText(f"{ideal_di_value:.2f}")
        
        QMessageBox.information(self, "Result", 
                               f"Calculated critical diameter:\n"
                               f"Actual: {di_value:.2f} mm\n"
                               f"Ideal: {ideal_di_value:.2f} mm")
        
        # Update plot
        self.update_jominy_plot()
        
        # Update comparison plot with DI from Jominy curve
        self.update_di_comparison()
    
    def update_jominy_plot(self):
        composition = self.get_composition()
        if composition is None:
            return
        
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, grain_size, hardness_after_quenching = composition
        
        # Calculate Just curve
        distances_just = np.linspace(0, 50, 100)  # Limit to 50 mm
        hardness_just = just_jominy_curve(C, Mn, Si, Cr, Ni, Mo, V, grain_size, distances_just)
        
        # Calculate ASTM curve
        hardness_astm = astm_jominy_curve(C, Mn, Si, Ni, Cr, Mo, Cu, grain_size, distances_just)
        
        # Get experimental data
        experimental_curves = self.get_jominy_data()
        
        # Update plot
        steel_name = self.steel_combo.currentText()
        
        # Pass all curves at once
        self.jominy_plot.plot_jominy(distances_just, hardness_just, distances_just, hardness_astm, 
                                    experimental_curves, steel_name, hardness_after_quenching, self.show_just_method)
    
    def save_composition(self):
        composition = self.get_composition()
        if composition is None:
            return
        
        C, Mn, Si, Cr, Ni, Mo, Cu, V, W, grain_size, hardness_after_quenching = composition
        
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Steel Composition", "", "JSON Files (*.json);;CSV Files (*.csv)")
        
        if file_path:
            if file_path.endswith('.json'):
                data = {
                    "C": C, "Mn": Mn, "Si": Si, "Cr": Cr, "Ni": Ni, 
                    "Mo": Mo, "Cu": Cu, "V": V, "W": W, 
                    "grain_size": grain_size, "hardness_after_quenching": hardness_after_quenching
                }
                with open(file_path, 'w') as f:
                    json.dump(data, f)
            elif file_path.endswith('.csv'):
                with open(file_path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(["Component", "Value [%]"])
                    writer.writerow(["C", C])
                    writer.writerow(["Mn", Mn])
                    writer.writerow(["Si", Si])
                    writer.writerow(["Cr", Cr])
                    writer.writerow(["Ni", Ni])
                    writer.writerow(["Mo", Mo])
                    writer.writerow(["Cu", Cu])
                    writer.writerow(["V", V])
                    writer.writerow(["W", W])
                    writer.writerow(["grain_size", grain_size])
                    writer.writerow(["hardness_after_quenching", hardness_after_quenching])
    
    def load_composition(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Steel Composition", "", "JSON Files (*.json);;CSV Files (*.csv)")
        
        if file_path:
            if file_path.endswith('.json'):
                with open(file_path, 'r') as f:
                    data = json.load(f)
                    self.c_input.setText(str(data.get("C", 0)))
                    self.mn_input.setText(str(data.get("Mn", 0)))
                    self.si_input.setText(str(data.get("Si", 0)))
                    self.cr_input.setText(str(data.get("Cr", 0)))
                    self.ni_input.setText(str(data.get("Ni", 0)))
                    self.mo_input.setText(str(data.get("Mo", 0)))
                    self.cu_input.setText(str(data.get("Cu", 0)))
                    self.v_input.setText(str(data.get("V", 0)))
                    self.w_input.setText(str(data.get("W", 0)))
                    self.grain_size_input.setText(str(data.get("grain_size", 7)))
                    self.hardness_after_quenching_input.setText(str(data.get("hardness_after_quenching", 50)))
            elif file_path.endswith('.csv'):
                with open(file_path, 'r') as f:
                    reader = csv.reader(f)
                    next(reader)  # Skip header
                    for row in reader:
                        if len(row) >= 2:
                            component = row[0].strip()
                            value = row[1].strip()
                            if component == "C":
                                self.c_input.setText(value)
                            elif component == "Mn":
                                self.mn_input.setText(value)
                            elif component == "Si":
                                self.si_input.setText(value)
                            elif component == "Cr":
                                self.cr_input.setText(value)
                            elif component == "Ni":
                                self.ni_input.setText(value)
                            elif component == "Mo":
                                self.mo_input.setText(value)
                            elif component == "Cu":
                                self.cu_input.setText(value)
                            elif component == "V":
                                self.v_input.setText(value)
                            elif component == "W":
                                self.w_input.setText(value)
                            elif component == "grain_size":
                                self.grain_size_input.setText(value)
                            elif component == "hardness_after_quenching":
                                self.hardness_after_quenching_input.setText(value)
        
        # Update plot after loading composition
        self.update_jominy_plot()
    
    def save_jominy_data(self):
        experimental_curves = self.get_jominy_data()
        
        if not experimental_curves:
            QMessageBox.warning(self, "Error", "No data to save.")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Jominy Data", "", "CSV Files (*.csv)")
        
        if file_path:
            with open(file_path, 'w', newline='') as f:
                writer = csv.writer(f)
                headers = ["Distance [mm]"]
                for i in range(4):
                    if self.curve_checkboxes[i].isChecked():
                        headers.append(f"Hardness {i+1} [HRC]")
                writer.writerow(headers)
                
                for row in range(self.jominy_table.rowCount()):
                    dist_item = self.jominy_table.item(row, 0)
                    if dist_item and dist_item.text():
                        row_data = [dist_item.text()]
                        for i in range(4):
                            if self.curve_checkboxes[i].isChecked():
                                hard_item = self.jominy_table.item(row, i+1)
                                row_data.append(hard_item.text() if hard_item and hard_item.text() else "")
                        writer.writerow(row_data)
    
    def load_jominy_data(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Jominy Data", "", "CSV Files (*.csv)")
        
        if file_path:
            # Clear table
            self.jominy_table.setRowCount(0)
            
            with open(file_path, 'r') as f:
                reader = csv.reader(f)
                headers = next(reader)  # Read headers
                
                # Determine which columns to load
                load_columns = [0]  # Always load distance
                for i, header in enumerate(headers[1:], 1):
                    if "Hardness" in header:
                        load_columns.append(i)
                
                # Load data
                for row in reader:
                    if len(row) >= len(headers):
                        row_index = self.jominy_table.rowCount()
                        self.jominy_table.insertRow(row_index)
                        
                        for col_idx, data_col in enumerate(load_columns):
                            if data_col < len(row):
                                self.jominy_table.setItem(row_index, col_idx, QTableWidgetItem(row[data_col]))
            
            # Update plot
            self.update_jominy_plot()
    
    def save_jominy_plot(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Jominy Plot", "", "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg)")
        
        if file_path:
            self.jominy_plot.fig.savefig(file_path)
    
    def save_results(self):
        if not self.di_results:
            QMessageBox.warning(self, "Error", "No results to save. Perform calculations first.")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Results", "", "CSV Files (*.csv);;Text Files (*.txt)")
        
        if file_path:
            with open(file_path, 'w', newline='') as f:
                if file_path.endswith('.csv'):
                    writer = csv.writer(f)
                    writer.writerow(["Method", "DI [mm]"])
                    for method, value in self.di_results.items():
                        if method != "Steel":
                            writer.writerow([method, f"{value:.2f}"])
                    
                    if self.jominy_di is not None:
                        writer.writerow(["Experimental (actual)", f"{self.jominy_di:.2f}"])
                        writer.writerow(["Experimental (ideal)", f"{self.ideal_di:.2f}"])
                    
                    writer.writerow(["Steel", self.di_results.get("Steel", "No data")])
                    writer.writerow(["Calculation date", datetime.now().strftime("%Y-%m-%d %H:%M:%S")])
                else:
                    f.write("CRITICAL DIAMETER DI CALCULATION RESULTS\n")
                    f.write("========================================\n\n")
                    
                    for method, value in self.di_results.items():
                        if method != "Steel":
                            f.write(f"{method}: {value:.2f} mm\n")
                    
                    if self.jominy_di is not None:
                        f.write(f"\nExperimental (actual): {self.jominy_di:.2f} mm\n")
                        f.write(f"Experimental (ideal): {self.ideal_di:.2f} mm\n")
                    
                    f.write(f"\nSteel: {self.di_results.get('Steel', 'No data')}\n")
                    f.write(f"Calculation date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SteelCalculatorApp()
    window.show()
    sys.exit(app.exec_())
