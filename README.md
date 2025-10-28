Steel Heat Treatment Analysis Suite
A comprehensive collection of professional metallurgical software tools for steel heat treatment analysis, prediction, and simulation. This suite includes three specialized applications for tempering analysis, transformation diagrams, and hardenability calculations.

📋 Applications Overview
1. Temperdata 1.0
Steel Tempering Analysis Software

Temperdata is designed for visualization and analysis of steel tempering process data. It enables viewing tempering curves for various steel grades, reviewing chemical compositions, and exporting data and plots.

Key Features:
Data Visualization: Interactive tempering curves with tooltips showing detailed information

Steel Database: Comprehensive database of steel grades with chemical compositions

Multi-tab Interface: Plot visualization, steel information, data tables, and user manual

Data Export: Save plots as images and data as CSV files

Auto-load Functionality: Automatic loading of default data files

Interactive Hover: Detailed point information on plot hover

Supported Data:
Tempering temperature ranges

Multiple tempering times

Final hardness measurements

Complete chemical compositions

Data source references

2. Quenching Studio 1.0
Steel Transformation Analysis Software

Quenching Studio provides advanced analysis of phase transformations in steels during heat treatment. Based on scientific models, it calculates critical temperatures, transformation kinetics, and predicts final properties.

Key Features:
TTT Diagrams: Time-Temperature-Transformation diagrams

CCT Diagrams: Continuous-Cooling-Transformation diagrams

Phase Transformation Analysis: Ferritic, pearlitic, bainitic, and martensitic transformations

Hardness Prediction: Vickers hardness as function of cooling rate

Steel Database: Pre-loaded with common steel grades

Multiple Calculation Methods: ASTM/Grossman, Just, de Cremona methods

Professional Visualization: High-quality plots and comparison charts

Technical Capabilities:
Critical temperature calculations (Ae1, Ae3, Bs, Ms)

Phase fraction evolution during cooling

Hardness prediction based on cooling rate

Support for custom steel compositions

Grain size effects on transformations

3. Jomina Analyzer 1.0
Professional Steel Hardenability Analysis Tool

Jomina Analyzer specializes in Jominy end-quench test analysis and hardenability predictions. It provides multiple calculation methods for critical diameter and Jominy curve prediction.

Key Features:
Multiple DI Methods: ASTM/Grossman, Just (1969), de Cremona (1970)

Jominy Curve Prediction: Calculate hardenability curves from composition

Experimental Data Support: Up to 4 experimental curves with comparison

Critical Diameter Calculation: Actual and ideal DI values

Steel Database: Comprehensive material database

Professional Interface: Multi-tab layout with advanced visualization

Calculation Methods:
ASTM/Grossman method for critical diameter

Just method for DI and Jominy curves

de Cremona method with comprehensive alloy factors

ASTM A255 standard compliance

Half-martensitic hardness calculations

🚀 Installation
Prerequisites
Python 3.7 or higher

Required packages:

text
PyQt5
pandas
numpy
matplotlib
scipy
Installation Steps
Clone the repository:


Install required packages:

bash
pip install PyQt5 pandas numpy matplotlib scipy
Run individual applications:

bash
python "Temperdata 1.0 EN.py"
python "Quenching studio 1.0 EN.py"
python "Jomina Analyser EN 1.0.py"
📖 Usage
Temperdata 1.0
Launch the application

Load CSV data (automatic or manual)

Select steel grade from dropdown

Navigate through tabs:

Plot: View tempering curves

Steel Information: Chemical composition

Data Table: Tabular data view

Manual: User documentation

Quenching Studio 1.0
Select steel from database or input custom composition

Set temperature parameters and cooling rates

Choose calculation type:

Single cooling rate analysis

Hardness vs cooling rate

TTT diagram generation

CCT diagram generation

View results in respective tabs

Jomina Analyzer 1.0
Select steel grade or input composition

Calculate critical diameter using multiple methods

Input experimental Jominy data (optional)

Compare calculated vs experimental curves

Analyze DI method comparisons

🔬 Technical Details
Data Sources
Tempering Data: Based on Kaggle database of tempering data for carbon and low alloy steels

Transformation Models: Scientific models from published research

Steel Database: Comprehensive collection of common industrial steels

Calculation Methods
ASTM/Grossman: Industry standard for hardenability

Just Method: German methodology (1969)

de Cremona: French approach (1970)

Empirical Models: Based on Maynier et al. hardness equations

Supported Steel Types
Carbon steels

Low-alloy steels

Case-hardening steels

Heat-treatable steels

Nitriding steels

📊 Output Features
Graphical Outputs
Professional-grade matplotlib plots

Interactive tooltips and hover information

Multiple curve comparisons

Export to PNG, PDF, SVG formats

Data Export
CSV format for numerical data

JSON for composition data

High-resolution images

Comprehensive results tables

🔧 System Requirements
OS: Windows, Linux, macOS

Memory: 4GB RAM minimum

Display: 1280x720 resolution minimum (Full HD recommended)

Storage: 100MB free space

📜 License
This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

Permissions
✅ Commercial use

✅ Modification

✅ Distribution

✅ Patent use

✅ Private use

Conditions
📝 License and copyright notice

📝 State changes

📝 Disclose source

📝 Same license

Limitations
❌ Liability

❌ Warranty

🤝 Contributing
We welcome contributions from the metallurgical and materials science community:

Fork the repository

Create a feature branch

Commit your changes

Push to the branch

Create a Pull Request

Areas for Contribution
Additional steel grades

Improved calculation models

User interface enhancements

Documentation improvements

Bug fixes and optimizations

📞 Support and Contact
Primary Contact: Marek Góral
Email: m_goral@interia.pl

Technical Support
Create an issue on GitHub

Email for direct support

Documentation available in each application

🎓 Educational Use
This software suite is particularly valuable for:

Universities: Materials science and engineering courses

Research Institutions: Heat treatment studies

Industry Professionals: Process optimization and analysis

Students: Learning steel transformation behavior

🔍 Scientific Background
The applications are based on established metallurgical principles and empirical models from:

ASTM A255 standard

Grange and Andrews equations for critical temperatures

Van Bohemen martensite transformation models

Maynier hardness equations

Scientific publications in Materials Science and Technology



🙏 Acknowledgments
Original Codebase: Transformation diagrams code from arthursn/transformation-diagrams

Data Sources: Kaggle tempering data database

Testing: Rzeszow University of Technology

Development: Regional Excellence Initiative program

Disclaimer: This software is intended for educational and professional analysis purposes. While based on established scientific models, actual heat treatment results may vary based on specific conditions and material variations. Always verify critical calculations with experimental data.
