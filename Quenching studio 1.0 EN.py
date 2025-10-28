import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from abc import abstractmethod
from scipy import integrate
from scipy.interpolate import splrep, splev, interp1d
from scipy.optimize import root
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QGridLayout, QLabel, QLineEdit, 
                             QPushButton, QGroupBox, QDoubleSpinBox, QTextEdit,
                             QTabWidget, QComboBox, QScrollArea, QMessageBox,
                             QTableWidget, QTableWidgetItem, QHeaderView,
                             QSplitter, QFrame, QSizePolicy, QFileDialog)
from PyQt5.QtCore import Qt, QSize
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import cm
import json

"""
Quenching Studio - Steel Transformation Analysis Software

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

This software is based on code from the repository:
https://github.com/arthursn/transformation-diagrams/tree/master

Contact: m_goral@interia.pl
"""

R = 8.314459
K = 273.15

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
    """
    Convert temperature from Fahrenheit to Celsius
    """
    return (TF - 32.)*5./9.

def FC(**comp):
    """
    Function expressing the influence of alloying elements on ferritic transformation kinetics
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return np.exp((1.0 + 6.31*C + 1.78*Mn + 0.31*Si + 1.12*Ni + 2.7*Cr + 4.06*Mo))

def PC(**comp):
    """
    Function expressing the influence of alloying elements on pearlitic transformation kinetics
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return np.exp(-4.25 + 4.12*C + 4.36*Mn + 0.44*Si + 1.71*Ni + 3.33*Cr + 5.19*np.sqrt(Mo))

def BC(**comp):
    """
    Function expressing the influence of alloying elements on bainitic transformation kinetics
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return np.exp(-10.23 + 10.18*C + 0.85*Mn + 0.55*Ni + 0.9*Cr + 0.36*Mo)

def Ae1_Grange(**comp):
    """
    Grange equation for Ae1
    """
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    return FahrenheitToCelsius(1333 - 25*Mn + 40*Si - 26*Ni + 42*Cr)

def Ae3_Grange(**comp):
    """
    Grange equation for Ae3
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    return FahrenheitToCelsius(1570 - 323*C - 25*Mn + 80*Si - 32*Ni - 3*Cr)

def Ae1_Andrews(**comp):
    """
    Andrews equation for Ae1
    """
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    W = comp.get('W', 0)
    As = comp.get('As', 0)
    return 723 - 16.9*Ni + 29.1*Si + 6.38*W - 10.7*Mn + 16.9*Cr + 290*As

def Ae3_Andrews(**comp):
    """
    Andrews equation for Ae3
    """
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
    """
    Bainite start temperature from Li's work
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 637 - 58*C - 35*Mn - 15*Ni - 34*Cr - 41*Mo

def Bs_VanBohemen(**comp):
    """
    Bainite start temperature from Van Bohemen's work
    [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 839 - (86*Mn + 23*Si + 67*Cr + 33*Ni + 75*Mo) - 270*(1 - np.exp(-1.33*C))

def Ms_Andrews(**comp):
    """
    Andrews equation for Ms
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    Co = comp.get('Co', 0)
    return 539 - 423*C - 30.4*Mn - 17.7*Ni - 12.1*Cr - 7.5*Mo + 10*Co - 7.5*Si

def alpha_martensite_VanBohemen(**comp):
    """
    Martensitic transformation rate constant
    [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 1e-3*(27.2 - (0.14*Mn + 0.21*Si + 0.11*Cr + 0.08*Ni + 0.05*Mo) - 19.8*(1-np.exp(-1.56*C)))

def Ms_VanBohemen(**comp):
    """
    Martensite start temperature
    [1] S.M.C. van Bohemen, Mater. Sci. Technol. 28 (2012) 487–495.
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return 565 - (31*Mn + 13*Si + 10*Cr + 18*Ni + 12*Mo) - 600*(1-np.exp(-0.96*C))

def Hv_martensite(phi700, **comp):
    """
    Empirical Vickers hardness equation for martensite
    (Maynier et al.)
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    return 127 + 949*C + 27*Si + 11*Mn + 8*Ni + 16*Cr + 21*np.log10(phi700*3600)

def Hv_bainite(phi700, **comp):
    """
    Empirical Vickers hardness equation for bainite
    (Maynier et al.)
    """
    C = comp.get('C', 0)
    Mn = comp.get('Mn', 0)
    Si = comp.get('Si', 0)
    Ni = comp.get('Ni', 0)
    Cr = comp.get('Cr', 0)
    Mo = comp.get('Mo', 0)
    return -323 + 185*C + 330*Si + 153*Mn + 65*Ni + 144*Cr + 191*Mo + \
        (89 + 53*C - 55*Si - 22*Mn - 10*Ni - 20*Cr - 33*Mo)*np.log10(phi700*3600)

def Hv_ferrite_pearlite(phi700, **comp):
    """
    Empirical Vickers hardness equation for ferrite + pearlite
    (Maynier et al.)
    """
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
    """
    Alloy properties (composition in wt.% and austenite grain size)
    """

    def __init__(self, gs, **w):
        # Grain size
        self.gs = gs

        # Alloy composition
        self.w = w
        # Main elements
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

        # Hardness
        self.Hv_martensite = lambda phi700: Hv_martensite(phi700, **w)
        self.Hv_bainite = lambda phi700: Hv_bainite(phi700, **w)
        self.Hv_ferrite_pearlite = lambda phi700: Hv_ferrite_pearlite(phi700, **w)

    def format_composition(self, vmin=0):
        fmt = []
        for k, v in self.w.items():
            if v > vmin:
                fmt.append('{:g}{:}'.format(v, k))
        fmt.insert(0, 'Fe')  # assumes it's steel
        fmt = '-'.join(fmt) + ' (wt.%)'
        fmt += '\nASTM grain size {:}'.format(self.gs)
        return fmt

class SigmoidalFunction(object):
    """
    Abstract class for S(X) and I(X) functions. After initialization,
    computes function values for a given range [xmin, xmax]
    and creates a spline interpolator. Returned values are computed by
    the interpolator. This method has the advantage of being able to process x as an array
    (or any other kind of iterator)
    """
    tck = None  # knots, coefficients and spline degree
    tck_inv = None  # inverse function spline parameters

    def __new__(cls, x):
        """
        __new__ behavior changed to return interpolated
        function value
        """
        if cls is SigmoidalFunction:
            raise TypeError("Cannot instantiate abstract class SigmoidalFunction")

        # Check required subclass attributes
        for var in ['xmin', 'xmax', 'ymin', 'ymax', 'n']:
            if not hasattr(cls, var):
                raise NotImplementedError(
                    'Class {} does not have required attribute `{}`'.format(cls, var))

        # Here S(X) or I(X) is returned
        return cls.val(x)

    @staticmethod
    def f(x):
        """
        Function to integrate
        """
        pass

    @classmethod
    def val(cls, x):
        """
        Computes SigmoidalFunction(x)
        """
        if hasattr(x, '__iter__') and not isinstance(x, str):
            x = np.array(x)
            xmin = x[x > 0].min()
        else:
            xmin = x
        xmin = min(cls.xmin, xmin)

        # initialize spline if not yet initialized or if xmin is lower than lower bound
        if xmin < cls.xmin or cls.tck is None:
            cls.xmin = xmin
            cls.init_spline()

        return splev(x, cls.tck)

    @classmethod
    def inv(cls, y):
        """
        Computes inverse function SigmoidalFunction^-1(y)
        """
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
        """
        Initializes spline
        """
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

    """
    S(X) function computed using numerical integration and
    spline interpolation
    """
    @staticmethod
    def f(x):
        return 1./(x**(0.4*(1. - x))*(1. - x)**(0.4*x))

class I(SigmoidalFunction):
    """
    I(X) function computed using numerical integration and
    spline interpolation
    """
    n = 999
    xmin = 0.001
    xmax = 0.999
    ymin = 0.29961765
    ymax = 4.05928646

    @staticmethod
    def f(x):
        return 1./(x**(2.*(1. - x)/3.)*(1. - x)**(2.*x/3.))

class PhaseTransformation(object):
    """
    Abstract class for calculating diffusion phase transformation kinetics
    """

    def __init__(self, alloy):
        self.alloy = alloy
        self.initialize()
        # Check required object attributes
        for var in ['comp_factor', 'Ts', 'Tf', 'Hv']:
            if not hasattr(self, var):
                raise NotImplementedError(
                    'Object {} does not have required attribute `{}`'.format(self, var))

    @classmethod
    def __init_subclass__(cls):
        # Check required subclass attributes
        for var in ['Q', 'n1', 'n2']:
            if not hasattr(cls, var):
                raise NotImplementedError(
                    'Class {} does not have required attribute `{}`'.format(cls, var))

    @abstractmethod
    def initialize(self):
        pass

    def get_transformation_factor(self, T):
        """
        Computes transformation factor for given temperature T

        Parameters
        ----------
        T : float or iterable
            Temperature. Can be provided as array

        Returns
        -------
        F : float or iterable
            Transformation factor with same shape as T
        """
        return self.comp_factor/(2**(self.n1*self.alloy.gs)*(self.Ts - T)**self.n2*np.exp(-self.Q/(R*(T + K))))

    def get_transformation_time(self, T, f):
        """
        Computes time required to transform material to fraction f at temperature T

        Parameters
        ----------
        T : float or iterable
            Temperature. Can be provided as array
        f : float or iterable
            Transformed fraction. If iterable, must have same shape as T

        Returns
        -------
        t : float or iterable
            Transformation time with same shape as T
        """
        return S(f)*self.get_transformation_factor(T)

    def get_transformation_temperature(self, Tini, Tfin, cooling_rate, f, dT=1.0):
        """
        Computes temperature at which material transforms to fraction f
        during cooling from Tini to Tfin with cooling rate

        Parameters
        ----------
        Tini : float
            Initial temperature
        Tfin : float
            Final temperature
        cooling_rate : float or iterable
            Cooling rate(s). Can be provided as array
        f : float or iterable
            Transformed fraction. If iterable, must have same shape as cooling_rate

        Returns
        -------
        T : float or iterable
            Transformation temperature with same shape as cooling_rate
        """
        dt = dT/np.array(cooling_rate)
        nt = len(dt) if hasattr(dt, '__len__') else 1
        T = np.arange(Tini, Tfin, -dT)
        nucleation_time = np.full((nt, len(T)), 0, dtype=float)

        filtr = T < self.Ts
        nucleation_time[:, filtr] = np.outer(dt, 1./self.get_transformation_factor(T[filtr]))
        nucleation_time = nucleation_time.cumsum(axis=1)

        Tt = np.full(nt, np.nan, dtype=float)  # Transformation temperature for given fraction f

        # Check indices where nucleation_time is greater than threshold S(f)
        # First occurrence is transformation temperature
        Sf = S(f)
        for i, n_time in enumerate(nucleation_time):
            idx, = np.where(n_time >= Sf)
            if len(idx) > 0:
                Tt[i] = T[idx[0]]

        return float(Tt) if nt == 1 else Tt

    def get_transformed_fraction(self, t, T, n=1000):
        """
        Computes transformed fraction for given thermal cycle T(t)

        Parameters
        ----------
        t : iterable
            Time
        T : iterable
            Temperatures at time instances t
        n : int (optional)
            Number of points where transformed fractions are computed
            Default: 1000

        Returns
        -------
        t, T, f : tuple
            Tuple with arrays of time, temperature, phase fraction computed
            at n points
        """
        if len(t) > 3:
            # Fit T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Use linear interpolator
            t2T = interp1d(t, T)

        # To ensure algorithm convergence, thermal cycle T(t) is
        # adjusted by spline, and nucleation time is computed by
        # increments dt = (max(t) - min(t))/n
        dt = (max(t) - min(t))/(n - 1)
        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        nucleation_time = np.full(t.shape, 0, dtype=float)
        f = np.full(T.shape, 0, dtype=float)

        # Compute nucleation time only for T lower than transformation
        # start temperature and higher than Tf
        filtr = (T < self.Ts) & (T > self.Tf)
        if np.any(filtr):
            nucleation_time[filtr] = dt/self.get_transformation_factor(T[filtr])
            nucleation_time = nucleation_time.cumsum()
            if T[0] < self.Ts:
                # This is the factor corresponding to transformed fraction at t[0]
                nucleation_time += min(t)/self.get_transformation_factor(T[0])

            # New filter: compute f only for nucleation_time within
            # S.inv(y) limits
            filtr = (nucleation_time >= S.ymin) & (nucleation_time <= S.ymax)
            if np.any(filtr):
                f[filtr] = S.inv(nucleation_time[filtr])
                f[nucleation_time < S.ymin] = 0
                f[nucleation_time > S.ymax] = 1

        return t, T, f

class Ferrite(PhaseTransformation):
    """
    Austenite to ferrite transformation
    """
    Q = 27500*4.184  # activation energy
    n1 = 0.41  # exponential coefficient 1
    n2 = 3  # exponential coefficient 2

    def initialize(self):
        self.comp_factor = self.alloy.FC  # composition factor for computing transformation time
        self.Ts = self.alloy.Ae3  # transformation start temperature
        self.Tf = self.alloy.Bs  # transformation finish temperature
        self.Hv = self.alloy.Hv_ferrite_pearlite

class Pearlite(PhaseTransformation):
    """
    Austenite to pearlite transformation
    """
    Q = 27500*4.184  # activation energy
    n1 = 0.32  # exponential coefficient 1
    n2 = 3  # exponential coefficient 2

    def initialize(self):
        self.comp_factor = self.alloy.PC  # composition factor for computing transformation time
        self.Ts = self.alloy.Ae1  # transformation start temperature
        self.Tf = self.alloy.Bs  # transformation finish temperature
        self.Hv = self.alloy.Hv_ferrite_pearlite

class Bainite(PhaseTransformation):
    """
    Austenite to bainite transformation
    """
    Q = 27500*4.184  # activation energy
    n1 = 0.29  # exponential coefficient 1
    n2 = 2  # exponential coefficient 2

    def initialize(self):
        self.comp_factor = self.alloy.BC  # composition factor for computing transformation time
        self.Ts = self.alloy.Bs  # transformation start temperature
        self.Tf = self.alloy.Ms  # transformation finish temperature
        self.Hv = self.alloy.Hv_bainite

class Martensite:
    """
    Diffusionless austenite to martensite transformation
    """

    def __init__(self, alloy):
        self.alloy = alloy
        self.Ts = self.alloy.Ms
        self.Hv = self.alloy.Hv_martensite

    def get_transformed_fraction(self, t, T, n=1000):
        """
        Computes transformed martensite fraction for given thermal cycle T(t)
        using Koistinen-Marburger equation

        Parameters
        ----------
        t : iterable
            Time
        T : iterable
            Temperatures at time instances t
        n : int (optional)
            Number of points where transformed fractions are computed
            Default: 1000

        Returns
        -------
        t, T, f : tuple
            Tuple with arrays of time, temperature, phase fraction computed
            at n points
        """
        if len(t) > 3:
            # Fit T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Use linear interpolator
            t2T = interp1d(t, T)

        t = np.linspace(min(t), max(t), n)
        T = t2T(t)
        f = np.full(T.shape, 0, dtype=float)

        filtr = T < self.alloy.Ms
        if np.any(filtr):
            f[filtr] = 1 - np.exp(-self.alloy.alpha_martensite*(self.alloy.Ms - T[filtr]))
        return t, T, f

class TransformationDiagrams:
    """
    Transformation diagrams class
    """

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
        """
        Computes transformation curves for given thermal cycle T(t)

        Parameters
        ----------
        t : iterable
            Time
        T : iterable
            Temperatures at time instances t
        n : int (optional)
            Number of points where transformed fractions are computed
            Default: 1000

        Returns
        -------
        f : pandas DataFrame
            DataFrame containing time, temperature and phase fractions
            of ferrite, pearlite, bainite, martensite and austenite at n
            points, and Vickers hardness for each data point
        """
        # Uncorrected phase fractions
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

        # Create DataFrame with appropriate data type to avoid warnings
        f = pd.DataFrame({
            't': t.astype(float),
            'T': T.astype(float),
            'ferrite': np.zeros(len(t), dtype=float),
            'pearlite': np.zeros(len(t), dtype=float),
            'bainite': np.zeros(len(t), dtype=float),
            'martensite': np.zeros(len(t), dtype=float),
            'austenite': np.zeros(len(t), dtype=float)
        })

        # Set initial values
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
            x0 = [
                f.loc[i-1, 'ferrite'], 
                f.loc[i-1, 'pearlite'],
                f.loc[i-1, 'bainite'], 
                f.loc[i-1, 'martensite']
            ]

            # Solve nonlinear equation system to get corrected phase fractions
            res = root(lambda x: [f1(i, *x), f2(i, *x), f3(i, *x), f4(i, *x)], x0=x0)

            # Make sure we assign float values to avoid type warnings
            f.loc[i, 'ferrite'] = float(res.x[0])
            f.loc[i, 'pearlite'] = float(res.x[1])
            f.loc[i, 'bainite'] = float(res.x[2])
            f.loc[i, 'martensite'] = float(res.x[3])
            f.loc[i, 'austenite'] = float(1.0 - sum(res.x))

        phi700 = None

        try:
            T2t = interp1d(T, t)
            # Get cooling rate at 700 °C
            phi700 = 2./(T2t(699.) - T2t(701.))
            if phi700 == 0:
                phi700 = None
        except ValueError:
            # May occur for isothermal heat treatments
            pass

        if phi700 is not None:
            f['Hv'] = f['martensite']*self.martensite.Hv(phi700) + f['bainite']*self.bainite.Hv(phi700) + \
                (f['ferrite'] + f['pearlite'])*self.ferrite.Hv(phi700)
        else:
            f['Hv'] = np.nan

        return f.round(12)

    def draw_thermal_cycle(self, ax, t, T, n=100, **kwargs):
        """
        Draws thermal cycle (cooling curve) on AxesSubplot object

        Parameters
        ----------
        ax : AxesSubplot object
            Axis on which to draw thermal cycle curve
        t : iterable
            Time
        T : iterable
            Temperatures at time instances t
        n : int (optional)
            Number of points where T(t) points are evaluated by
            interpolation
            Default: 100

        Returns
        -------
        line : Line2D object
            Line2D object corresponding to drawn curve
        """

        if len(t) > 3:
            # Fit T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Use linear interpolator
            t2T = interp1d(t, T)

        t = np.linspace(min(t), max(t), n)
        T = t2T(t)

        kw = dict(color='k', ls='--')
        kw.update(kwargs)

        return ax.plot(t, T, **kw)

    def TTT(self, fs=1e-2, ff=.99, ax=None, **kwargs):
        """
        Draws TTT diagram

        Parameters
        ----------
        fs : float (optional)
            Phase fraction at transformation start
            Default: 1e-2 (1%)
        fs : float (optional)
            Phase fraction at transformation finish
            Default: .99 (99%)
        ax : AxesSubplot object (optional)
            Axis on which to draw TTT curve. If None, new axis is created
            Default: None
        **kwargs :
            Optional arguments passed to ax.plot(*args, **kwargs)

        Returns
        -------
        ax : AxesSubplot object
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()

        # Ferrite
        T = np.arange(self.alloy.Bs, self.alloy.Ae3)
        ts = self.ferrite.get_transformation_time(T, fs)  # start
        tf = self.ferrite.get_transformation_time(T, ff)  # finish
        ax.plot(ts, T, color=self.colors_dict['ferrite'],
                label='Ferrite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['ferrite'], ls='--',
                label='Ferrite {:.0f}%'.format(100*ff), **kwargs)
        df_ferrite = pd.DataFrame(dict(T_ferrite=T, ts_ferrite=ts, tf_ferrite=tf))

        # Pearlite
        T = np.arange(self.alloy.Bs, self.alloy.Ae1)
        ts = self.pearlite.get_transformation_time(T, fs)
        tf = self.pearlite.get_transformation_time(T, ff)
        ax.plot(ts, T, color=self.colors_dict['pearlite'],
                label='Pearlite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['pearlite'], ls='--',
                label='Pearlite {:.0f}%'.format(100*ff), **kwargs)
        df_pearlite = pd.DataFrame(dict(T_pearlite=T, ts_pearlite=ts, tf_pearlite=tf))

        # Bainite
        T = np.arange(self.alloy.Ms, self.alloy.Bs)
        ts = self.bainite.get_transformation_time(T, fs)
        tf = self.bainite.get_transformation_time(T, ff)
        ax.plot(ts, T, color=self.colors_dict['bainite'],
                label='Bainite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(tf, T, color=self.colors_dict['bainite'], ls='--',
                label='Bainite {:.0f}%'.format(100*ff), **kwargs)
        df_bainite = pd.DataFrame(dict(T_bainite=T, ts_bainite=ts, tf_bainite=tf))

        self.df_TTT = pd.concat([df_ferrite, df_pearlite, df_bainite], axis=1)

        # Draw Ae1 and Ae3 lines
        ax.axhline(self.alloy.Ae3, xmax=.1, color=self.colors_dict['ferrite'], ls=':')
        ax.axhline(self.alloy.Ae1, xmax=.1, color=self.colors_dict['pearlite'], ls=':')

        # Draw Bs and Ms lines
        ax.axhline(self.alloy.Bs, color=self.colors_dict['bainite'], ls=':')
        ax.axhline(self.alloy.Ms, color=self.colors_dict['martensite'])

        ax.set_xscale('log')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Temperature (°C)')
        ax.set_title(self.alloy.format_composition())

        xmin = ax.get_xlim()[0]
        ax.text(xmin*1.5, self.alloy.Ae3, 'Ae3',
                color=self.colors_dict['ferrite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ae1, 'Ae1',
                color=self.colors_dict['pearlite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Bs, 'Bs',
                color=self.colors_dict['bainite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ms, 'Ms',
                color=self.colors_dict['martensite'], ha='left', va='bottom')

        ax.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, -.15))
        fig.subplots_adjust(bottom=.2)

        return ax

    def CCT(self, Tini=900, fs=1e-2, ff=.99, phi_min=1e-4, phi_max=1e4, phi_steps=420, ax=None, **kwargs):
        """
        Draws CCT diagram

        Parameters
        ----------
        fs : float (optional)
            Phase fraction at transformation start
            Default: 1e-2 (1%)
        fs : float (optional)
            Phase fraction at transformation finish
            Default: .99 (99%)
        ax : AxesSubplot object (optional)
            Axis on which to draw TTT curve. If None, new axis is created
            Default: None
        **kwargs :
            Optional arguments passed to ax.plot(*args, **kwargs)

        Returns
        -------
        ax : AxesSubplot object
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()

        cooling_rates = 10**np.linspace(np.log10(phi_min), np.log10(phi_max), phi_steps)
        draw_cooling = kwargs.get('draw_cooling', True)

        # Ferrite
        Ts = self.ferrite.get_transformation_temperature(
            Tini, self.alloy.Bs, cooling_rates, fs)  # start
        Tf = self.ferrite.get_transformation_temperature(
            Tini, self.alloy.Bs, cooling_rates, ff)  # finish
        ax.plot(Ts/cooling_rates, Ts,
                color=self.colors_dict['ferrite'], label='Ferrite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['ferrite'],
                ls='--', label='Ferrite {:.0f}%'.format(100*ff), **kwargs)

        # Pearlite
        Ts = self.pearlite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, fs)
        Tf = self.pearlite.get_transformation_temperature(Tini, self.alloy.Bs, cooling_rates, ff)
        ax.plot(Ts/cooling_rates, Ts, color=self.colors_dict['pearlite'],
                label='Pearlite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['pearlite'],
                ls='--', label='Pearlite {:.0f}%'.format(100*ff), **kwargs)

        # Bainite
        Ts = self.bainite.get_transformation_temperature(Tini, self.alloy.Ms, cooling_rates, fs)
        Tf = self.bainite.get_transformation_temperature(Tini, self.alloy.Ms, cooling_rates, ff)
        ax.plot(Ts/cooling_rates, Ts,
                color=self.colors_dict['bainite'], label='Bainite {:.0f}%'.format(100*fs), **kwargs)
        ax.plot(Tf/cooling_rates, Tf, color=self.colors_dict['bainite'],
                ls='--', label='Bainite {:.0f}%'.format(100*ff), **kwargs)

        ax.axhline(self.alloy.Ae3, xmax=.1, color=self.colors_dict['ferrite'], ls=':')
        ax.axhline(self.alloy.Ae1, xmax=.1, color=self.colors_dict['pearlite'], ls=':')

        ax.axhline(self.alloy.Bs, color=self.colors_dict['bainite'], ls=':')
        ax.axhline(self.alloy.Ms, color=self.colors_dict['martensite'])

        # Draw cooling curves
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
        ax.text(xmin*1.5, self.alloy.Ae3, 'Ae3',
                color=self.colors_dict['ferrite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ae1, 'Ae1',
                color=self.colors_dict['pearlite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Bs, 'Bs',
                color=self.colors_dict['bainite'], ha='left', va='bottom')
        ax.text(xmin*1.5, self.alloy.Ms, 'Ms',
                color=self.colors_dict['martensite'], ha='left', va='bottom')

        ax.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, -.15))
        fig.subplots_adjust(bottom=.2)

        return ax

    def plot_phase_fraction(self, t, T, n=1000, xaxis='t', ax=None, **kwargs):
        """
        Draws phase fractions for given thermal cycle T(t)

        Parameters
        ----------
        t : iterable
            Time
        T : iterable
            Temperatures at time instances t
        n : int (optional)
            Number of points where T(t) points are evaluated by
            interpolation
            Default: 1000
        xaxis : string (optional)
            Variable to represent on x-axis. Possible values:
            'index', t', 'T', 'ferrite', 'pearlite', 'bainite', 'martensite',
            and 'austenite'
            Default: 't'
        ax : AxesSubplot object (optional)
            Axis on which to draw phase fraction curves. If None, then
            new axis is created
            Default: None

        Returns
        -------
        ax : AxesSubplot object
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        if len(t) > 3:
            # Fit T(t) by spline
            def t2T(t_): return splev(t_, splrep(t, T))
        else:
            # Use linear interpolator
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
            ax.plot(f[xaxis], f['martensite'],
                    color=self.colors_dict['martensite'], label='Martensite')
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
        self.setWindowTitle("Quenching studio 1.0")
        self.setGeometry(100, 50, 1800, 1000)  # Full HD size
        
        # Main widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
        # Create tab widget for all tabs
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)
        
        # Create input tab
        self.input_tab = QWidget()
        input_layout = QHBoxLayout(self.input_tab)  # Changed to QHBoxLayout for 3 columns
        
        # First column: Steel selection
        col1_layout = QVBoxLayout()
        steel_group = QGroupBox("Steel Selection")
        steel_layout = QVBoxLayout(steel_group)
        
        self.steel_combo = QComboBox()
        self.steel_combo.addItems(steel_database.keys())
        self.steel_combo.currentTextChanged.connect(self.load_steel_data)
        steel_layout.addWidget(QLabel("Select steel:"))
        steel_layout.addWidget(self.steel_combo)
        
        # Display steel information
        self.steel_info = QTextEdit()
        self.steel_info.setReadOnly(True)
        self.steel_info.setMaximumHeight(200)
        steel_layout.addWidget(QLabel("Steel information:"))
        steel_layout.addWidget(self.steel_info)
        
        col1_layout.addWidget(steel_group)
        col1_layout.addStretch()
        
        # Second column: Alloy composition
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
        
        # Save/Load composition buttons
        save_load_layout = QHBoxLayout()
        self.save_composition_btn = QPushButton("Save composition")
        self.save_composition_btn.clicked.connect(self.save_composition)
        self.load_composition_btn = QPushButton("Load composition")
        self.load_composition_btn.clicked.connect(self.load_composition)
        save_load_layout.addWidget(self.save_composition_btn)
        save_load_layout.addWidget(self.load_composition_btn)
        composition_layout.addLayout(save_load_layout, len(elements), 0, 1, 2)
        
        # Grain size input
        composition_layout.addWidget(QLabel("Grain Size"), len(elements)+1, 0)
        self.gs_input = QDoubleSpinBox()
        self.gs_input.setRange(1, 12)
        self.gs_input.setValue(7)
        self.gs_input.setDecimals(1)
        composition_layout.addWidget(self.gs_input, len(elements)+1, 1)
        
        col2_layout.addWidget(composition_group)
        col2_layout.addStretch()
        
        # Third column: Other elements
        col3_layout = QVBoxLayout()
        
        # Temperature parameters
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
        
        # Single cooling rate
        temp_layout.addWidget(QLabel("Single Cooling Rate (°C/s)"), 2, 0)
        self.single_cooling_input = QDoubleSpinBox()
        self.single_cooling_input.setRange(0.001, 10000)
        self.single_cooling_input.setValue(10)
        self.single_cooling_input.setDecimals(3)
        self.single_cooling_input.setSingleStep(1)
        temp_layout.addWidget(self.single_cooling_input, 2, 1)
        
        col3_layout.addWidget(temp_group)
        
        # Cooling rates for analysis
        cooling_group = QGroupBox("Cooling Rates for Analysis (°C/s)")
        cooling_layout = QVBoxLayout(cooling_group)
        
        self.cooling_rates_edit = QTextEdit()
        self.cooling_rates_edit.setPlainText(
            "1000\n300\n100\n30\n10\n3\n1\n0.3\n0.1\n0.03\n0.01\n0.003\n0.001"
        )
        self.cooling_rates_edit.setMaximumHeight(150)
        cooling_layout.addWidget(self.cooling_rates_edit)
        
        col3_layout.addWidget(cooling_group)
        
        # Calculation buttons
        self.calculate_single_btn = QPushButton("Calculate Single Cooling Rate")
        self.calculate_single_btn.clicked.connect(self.calculate_single_cooling)
        
        self.calculate_hardness_btn = QPushButton("Calculate Hardness vs Cooling Rate")
        self.calculate_hardness_btn.clicked.connect(self.calculate_hardness)
        
        self.calculate_ttt_btn = QPushButton("Calculate TTT Diagram")
        self.calculate_ttt_btn.clicked.connect(self.calculate_ttt)
        
        self.calculate_cct_btn = QPushButton("Calculate CCT Diagram")
        self.calculate_cct_btn.clicked.connect(self.calculate_cct)
        
        col3_layout.addWidget(self.calculate_single_btn)
        col3_layout.addWidget(self.calculate_hardness_btn)
        col3_layout.addWidget(self.calculate_ttt_btn)
        col3_layout.addWidget(self.calculate_cct_btn)
        col3_layout.addStretch()
        
        # Add columns to input layout
        input_layout.addLayout(col1_layout)
        input_layout.addLayout(col2_layout)
        input_layout.addLayout(col3_layout)
        
        # Add input tab to tab widget
        self.tab_widget.addTab(self.input_tab, "Input Parameters")
        
        # Single cooling rate tab
        self.single_tab = QWidget()
        single_layout = QVBoxLayout(self.single_tab)
        
        self.single_fig = Figure(figsize=(10, 8))
        self.single_canvas = FigureCanvas(self.single_fig)
        self.single_toolbar = NavigationToolbar(self.single_canvas, self)
        
        single_layout.addWidget(self.single_toolbar)
        single_layout.addWidget(self.single_canvas)
        
        # New tab: Single cooling rate table
        self.single_table_tab = QWidget()
        single_table_layout = QVBoxLayout(self.single_table_tab)
        
        self.single_table = QTableWidget()
        self.single_table.setColumnCount(2)
        self.single_table.setHorizontalHeaderLabels(["Property", "Value"])
        self.single_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        single_table_layout.addWidget(self.single_table)
        
        # Save to CSV button
        self.save_single_table_btn = QPushButton("Save results to CSV")
        self.save_single_table_btn.clicked.connect(self.save_single_table)
        single_table_layout.addWidget(self.save_single_table_btn)
        
        # Hardness vs cooling rate tab
        self.hardness_tab = QWidget()
        hardness_layout = QVBoxLayout(self.hardness_tab)
        
        self.fig1 = Figure(figsize=(10, 8))
        self.canvas1 = FigureCanvas(self.fig1)
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        
        hardness_layout.addWidget(self.toolbar1)
        hardness_layout.addWidget(self.canvas1)
        
        # New tab: Hardness vs cooling rate table
        self.hardness_table_tab = QWidget()
        hardness_table_layout = QVBoxLayout(self.hardness_table_tab)
        
        self.hardness_table = QTableWidget()
        self.hardness_table.setColumnCount(6)
        self.hardness_table.setHorizontalHeaderLabels(["Cooling Rate (°C/s)", "Ferrite", "Pearlite", "Bainite", "Martensite", "Hardness (HV)"])
        self.hardness_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        hardness_table_layout.addWidget(self.hardness_table)
        
        # Save to CSV button
        self.save_hardness_table_btn = QPushButton("Save results to CSV")
        self.save_hardness_table_btn.clicked.connect(self.save_hardness_table)
        hardness_table_layout.addWidget(self.save_hardness_table_btn)
        
        # TTT diagram tab
        self.ttt_tab = QWidget()
        ttt_layout = QVBoxLayout(self.ttt_tab)
        
        self.fig2 = Figure(figsize=(10, 8))
        self.canvas2 = FigureCanvas(self.fig2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        
        ttt_layout.addWidget(self.toolbar2)
        ttt_layout.addWidget(self.canvas2)
        
        # CCT diagram tab
        self.cct_tab = QWidget()
        cct_layout = QVBoxLayout(self.cct_tab)
        
        self.fig3 = Figure(figsize=(10, 8))
        self.canvas3 = FigureCanvas(self.fig3)
        self.toolbar3 = NavigationToolbar(self.canvas3, self)
        
        cct_layout.addWidget(self.toolbar3)
        cct_layout.addWidget(self.canvas3)
        
        # Instructions tab
        self.instruction_tab = QWidget()
        instruction_layout = QVBoxLayout(self.instruction_tab)
        
        instruction_text = QTextEdit()
        instruction_text.setReadOnly(True)
        instruction_text.setPlainText("""
This application is based on GNU code from the repository:<br>
        https://github.com/arthursn/transformation-diagrams/tree/master</p>

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
        
        # Add tabs to tab widget
        self.tab_widget.addTab(self.single_tab, "Single Cooling Rate")
        self.tab_widget.addTab(self.single_table_tab, "Single Cooling Rate - Table")
        self.tab_widget.addTab(self.hardness_tab, "Hardness vs Cooling Rate")
        self.tab_widget.addTab(self.hardness_table_tab, "Hardness vs Cooling Rate - Table")
        self.tab_widget.addTab(self.ttt_tab, "TTT Diagram")
        self.tab_widget.addTab(self.cct_tab, "CCT Diagram")
        self.tab_widget.addTab(self.instruction_tab, "Instructions")
        
        # Initialize variables
        self.alloy = None
        self.diagrams = None
        self.single_results = None
        self.hardness_results = None
        
        # Load initial steel data
        self.load_steel_data()
        
    def load_steel_data(self):
        steel_name = self.steel_combo.currentText()
        steel_data = steel_database[steel_name]
        
        # Update composition inputs
        for element, value in steel_data.items():
            if element in self.composition_inputs:
                self.composition_inputs[element].setValue(value)
        
        # Update grain size
        self.gs_input.setValue(steel_data.get("grain_size", 7))
        
        # Update steel information
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
                data = {
                    'composition': composition,
                    'grain_size': self.gs_input.value()
                }
                
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
                
                # Update composition inputs
                for element, value in composition.items():
                    if element in self.composition_inputs:
                        self.composition_inputs[element].setValue(value)
                
                # Update grain size
                self.gs_input.setValue(grain_size)
                
                # Set steel combobox to Custom
                self.steel_combo.setCurrentText("Custom")
                
                QMessageBox.information(self, "Success", "Composition loaded successfully!")
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load composition: {str(e)}")
        
    def calculate_single_cooling(self):
        try:
            # Get composition values
            composition = self.get_composition()
            
            # Get other parameters
            gs = self.gs_input.value()
            Tini = self.tini_input.value()
            Tfin = self.tfin_input.value()
            cooling_rate = self.single_cooling_input.value()
            
            # Create alloy and diagrams objects
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            
            # Calculate phase fractions and hardness
            total_time = (Tini - Tfin) / cooling_rate
            f = self.diagrams.get_transformed_fraction([0, total_time], [Tini, Tfin])
            f_fin = f.iloc[-1]
            
            # Save results
            self.single_results = {
                "Cooling Rate": cooling_rate,
                "Ferrite": f_fin['ferrite'],
                "Pearlite": f_fin['pearlite'],
                "Bainite": f_fin['bainite'],
                "Martensite": f_fin['martensite'],
                "Hardness": f_fin['Hv']
            }
            
            # Update single cooling rate plot
            self.update_single_cooling_plot(f, cooling_rate)
            
            # Update results table
            self.update_single_results_table(f_fin, cooling_rate)
            
            # Switch to single cooling rate tab
            self.tab_widget.setCurrentIndex(1)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
    
    def update_single_cooling_plot(self, f, cooling_rate):
        # Clear previous plots
        self.single_fig.clear()
        
        # Create subplots
        ax1 = self.single_fig.add_subplot(121)
        ax2 = self.single_fig.add_subplot(122)
        
        # Phase fractions vs time plot
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
        
        # Phase fractions pie chart
        phase_fractions = [f.iloc[-1][phase] for phase in phases]
        labels = [f'{phase.capitalize()}: {frac:.1%}' for phase, frac in zip(phases, phase_fractions) if frac > 0]
        sizes = [frac for frac in phase_fractions if frac > 0]
        colors = [color for phase, color in zip(phases, colors) if f.iloc[-1][phase] > 0]
        
        if sizes:
            ax2.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
            ax2.axis('equal')
            ax2.set_title('Final Phase Fractions')
        
        # Set main title
        self.single_fig.suptitle(self.alloy.format_composition())
        
        # Refresh canvas
        self.single_canvas.draw()
    
    def update_single_results_table(self, f_fin, cooling_rate):
        # Clear table
        self.single_table.setRowCount(6)
        self.single_table.setColumnCount(2)
        self.single_table.setHorizontalHeaderLabels(["Property", "Value"])
        
        # Set data
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
            # Get composition values
            composition = self.get_composition()
            
            # Get other parameters
            gs = self.gs_input.value()
            Tini = self.tini_input.value()
            Tfin = self.tfin_input.value()
            
            # Parse cooling rates
            cooling_text = self.cooling_rates_edit.toPlainText()
            cooling_rates = []
            for line in cooling_text.split('\n'):
                try:
                    cooling_rates.append(float(line.strip()))
                except ValueError:
                    pass
            
            # Create alloy and diagrams objects
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            
            # Calculate phase fractions and hardness
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
            
            # Save results
            self.hardness_results = pd.DataFrame({
                "Cooling Rate": cooling_rates,
                "Ferrite": f_ferr,
                "Pearlite": f_pear,
                "Bainite": f_bain,
                "Martensite": f_mart,
                "Hardness": Hv
            })
            
            # Update plots
            self.update_hardness_plots(cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv)
            
            # Update results table
            self.update_results_table(cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv)
            
            # Switch to hardness tab
            self.tab_widget.setCurrentIndex(3)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
    
    def update_hardness_plots(self, cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv):
        # Clear previous plots
        self.fig1.clear()
        
        # Create subplots
        ax1 = self.fig1.add_subplot(121)
        ax2 = self.fig1.add_subplot(122)
        
        # Phase fractions plot
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
        
        # Hardness plot
        ax2.plot(cooling_rates, Hv, 'b-')
        ax2.set_xlabel('Cooling Rate (°C/s)')
        ax2.set_ylabel('Vickers Hardness (HV)')
        ax2.set_xscale('log')
        ax2.set_title('Hardness vs Cooling Rate')
        ax2.grid(True)
        
        # Set main title
        self.fig1.suptitle(self.alloy.format_composition())
        
        # Refresh canvas
        self.canvas1.draw()
    
    def update_results_table(self, cooling_rates, f_ferr, f_pear, f_bain, f_mart, Hv):
        # Clear table
        self.hardness_table.setRowCount(len(cooling_rates))
        self.hardness_table.setColumnCount(6)
        self.hardness_table.setHorizontalHeaderLabels(["Cooling Rate (°C/s)", "Ferrite", "Pearlite", "Bainite", "Martensite", "Hardness (HV)"])
        
        # Set data
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
            # Get composition values
            composition = self.get_composition()
            
            # Get other parameters
            gs = self.gs_input.value()
            
            # Create alloy and diagrams objects
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            
            # Update TTT diagram
            self.update_ttt_diagram()
            
            # Switch to TTT tab
            self.tab_widget.setCurrentIndex(5)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
    
    def update_ttt_diagram(self):
        # Clear previous plots
        self.fig2.clear()
        
        # Create plot
        ax = self.fig2.add_subplot(111)
        
        # Draw TTT diagram
        self.diagrams.TTT(ax=ax)
        
        # Refresh canvas
        self.canvas2.draw()
    
    def calculate_cct(self):
        try:
            # Get composition values
            composition = self.get_composition()
            
            # Get other parameters
            gs = self.gs_input.value()
            Tini = self.tini_input.value()
            
            # Create alloy and diagrams objects
            self.alloy = Alloy(gs=gs, **composition)
            self.diagrams = TransformationDiagrams(self.alloy)
            
            # Update CCT diagram
            self.update_cct_diagram(Tini)
            
            # Switch to CCT tab
            self.tab_widget.setCurrentIndex(6)
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
    
    def update_cct_diagram(self, Tini):
        # Clear previous plots
        self.fig3.clear()
        
        # Create plot
        ax = self.fig3.add_subplot(111)
        
        # Draw CCT diagram
        self.diagrams.CCT(Tini=Tini, ax=ax)
        
        # Refresh canvas
        self.canvas3.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = SteelAnalysisApp()
    window.show()
    sys.exit(app.exec_())
