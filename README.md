# Guide to Using the Conservation Analyser
The Conservation Analyser is a Python-based graphical user interface (GUI) application designed for analysing sequence conservation from a given multiple sequence alignment (MSA) file. This guide covers setting up your environment, installing necessary tools, and using the application.
# Setting Up the Environment

1. Install Python

Ensure Python 3.x is installed on your system. Python can be downloaded from the official website. Follow the installation instructions for your operating system.

2. Create a Virtual Environment (Optional but Recommended)

A virtual environment is recommended to manage dependencies. Open your terminal or command prompt and navigate to your project directory:

#Navigate to your project directory

cd path/to/your/project

#Create a virtual environment named 'env'

python -m venv env

#Activate the virtual environment

#On Windows

env\Scripts\activate

#On macOS/Linux

source env/bin/activate

3. Install Required Python Packages

With your virtual environment activated, install the necessary Python packages using pip:

pip install PyQt5 matplotlib scipy

Installing MAFFT and Clustal Omega

# MAFFT

Windows: Download the installer from the MAFFT official website and follow the installation instructions.

macOS/Linux: MAFFT can typically be installed via package managers.

On macOS with Homebrew: brew install mafft

On Linux (Debian/Ubuntu): sudo apt-get install mafft

# Clustal Omega

Windows: Download from the Clustal Omega website and follow the setup instructions.

macOS/Linux: Install via package managers.

On macOS with Homebrew: brew install clustal-omega

On Linux (Debian/Ubuntu): sudo apt-get install clustalo

Using the Conservation Analyser

Start the Application: Navigate to the directory containing your ConservationAnalyzer.py script. Ensure your virtual environment is activated, and start the application by running: 

python3 ConservationAnalyser.py

Select Alignment File: Click the "Browse" button to open a file dialogue. Select your alignment file in FASTA format.

Set Parameters: Enter values for the conservation and mutation thresholds (0-1). Optionally, specify a reference sequence ID to focus the analysis.

Choose Aligner: Select either MAFFT or Clustal Omega as the aligner. Ensure the chosen aligner is installed and accessible from your system's PATH.

Analyse: Click the "Analyze" button to perform the sequence conservation analysis. The application will display conserved and mutated sites, overall conservation rate, and statistical significance (p-values).

View Results: The conserved and mutated sites will be listed in separate tables within the application. Additionally, a plot of conservation rates across the alignment will be displayed.

# Notes

The application assumes that input sequences are either pre-aligned or will be aligned using the chosen aligner.

The analysis and plotting functionalities require a graphical environment to display the results.

If running the application in an environment without a desktop GUI (e.g., a remote server), ensure you have access to an X11 forwarding setup or use the application on a local machine.

# Troubleshooting

Aligner Not Found: Ensure MAFFT and/or Clustal Omega is installed and correctly added to your system's PATH. You might need to restart your terminal or computer after installation.

Python or PyQt5 Errors: Make sure you have the correct versions of Python and PyQt5 installed. Reinstalling or updating these packages can resolve many issues.

Other Issues: Check the console output for error messages. They can provide valuable clues about what might be going wrong.

By following this guide, you should be able to set up your environment, run the Conservation Analyzer, and start analysing sequence conservation with ease.


