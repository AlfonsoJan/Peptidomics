# Peptidomics Spring Boot Project

The majority of the biological machinery is made up of proteins, which are just lengthy chains of connected amino acids. Although this relationship is exceedingly difficult to assess, in theory, the structure, dynamics, and function of such proteins follow from the individual building pieces.

An oligopeptide, which is made up of building units spanning multiple amino acids, is frequently preferable. Many similar oligopeptides are also known to serve significant purposes on their own. On the broader range of such peptides, not much is known. Their configurational qualities are one of them. The oligopeptide's conformational flexibility, the secondary structure annotations that correlate to those shapes, and their three-dimensional capabilities.

This software is the creation of something that essentially combines PCA and multilinear regression to translate the protein structure between its actual configurations and a low-dimensional abstract data space.

## Requirements
````
Java (min. version 8)
    Spring (version 3.0.2)
    Spring Security (version 2.0.4-RELEASE)

Python (min. version 3.9)
    Numpy (version 1.23.1)
    MDAnalysis (version 2.4.2)

JavaScript
    Toastr (version 2.1.4)
    JSmol (version 14.32)
    JQuery (version 3.6.3)
    Bulma Tooltip (version 1.2.0)
    Plotly (version 2.20.0)
````
## Installation
````
Clone this repository to your local machine:
    git clone https://github.com/AlfonsoJan/Peptidomics.git
````

Install the required Python dependencies:

    pip install -r requirements.txt

## Usage

For now open the project in IntelliJ

Set in the application.properties the python.path-name to your python3 shortcurt (eg, py, python, python3)

And run the main in the class: PeptidomicsWebAppApplication.

## Contributors

    John Busker (j.a.busker@st.hanze.nl)
    Wouter Zeevat (w.h.zeevat@st.hanze.nl)