README Summary
This project demonstrates a complete workflow for downloading, processing, and visualizing a protein structure with a point mutation. The key features include:

PDB Structure Download & Parsing:
Downloads the GLP1R structure (PDB ID: 6X1A) from the RCSB PDB database and parses it using BioPython’s PDBParser.

Sequence Extraction & Mutation:
Extracts the amino acid sequence from the first protein chain and creates a mutant sequence by replacing the amino acid at position 7 with Leucine (P7L mutation). Both wild-type and mutant sequences are saved as FASTA files.

3D Visualization:
Uses Matplotlib to create 3D visualizations of the protein backbone:

A plot for the wild-type structure.
A plot for the mutant structure with the mutation highlighted.
A side-by-side comparison figure for clear visual contrast.
Output:
All generated files (PDB, FASTA sequences, and images) are stored in the protein_data directory.

Dependencies:
Python 3, BioPython, NumPy, and Matplotlib are required to run the script.

This script serves as a template for protein structure analysis and visualization, useful for educational purposes or as a starting point for more advanced structural bioinformatics projects.