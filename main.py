import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import PPBuilder, protein_letters_3to1

# Create directory for files
if not os.path.exists('protein_data'):
    os.makedirs('protein_data')

# Function to download a PDB structure 
# We'll use GLP1R structure as example (PDB ID: 6X1A)
def get_structure(pdb_id, file_path):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir=file_path, file_format='pdb')
    return os.path.join(file_path, f"pdb{pdb_id.lower()}.ent")

# Download GLP1R structure
pdb_file = get_structure('6X1A', 'protein_data')

# Parse the PDB file
parser = PDBParser(QUIET=True)
structure = parser.get_structure('GLPR', pdb_file)

# Extract the sequence from the structure
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    sequence = pp.get_sequence()
    break  # Just take the first chain's sequence

# Define the wild-type and mutant sequences based on the extracted sequence
# For demonstration, we'll modify the actual structure's sequence
wt_sequence = sequence
print(f"Original sequence length: {len(wt_sequence)}")
print(f"First 20 amino acids: {wt_sequence[:20]}")

# Create the P7L mutation (assuming position 7 is accessible)
# For a real-world scenario, you'd need to ensure this is the correct position
mutant_sequence = list(str(wt_sequence))

# Display the original amino acid at position 7
original_aa = mutant_sequence[6]  # 0-indexed, so position 7 is index 6
print(f"Original amino acid at position 7: {original_aa}")

# Replace with Leucine (L)
mutant_sequence[6] = 'L'
mutant_sequence = ''.join(mutant_sequence)

# Save sequences to FASTA files
with open('protein_data/wt_glpr.fasta', 'w') as f:
    SeqIO.write(SeqRecord(Seq(str(wt_sequence)), id="wild_type_glpr", description=""), f, "fasta")

with open('protein_data/mutant_glpr.fasta', 'w') as f:
    SeqIO.write(SeqRecord(Seq(mutant_sequence), id="mutant_glpr_P7L", description=""), f, "fasta")

# Function to highlight specific residues
class MutationSelect(Select):
    def __init__(self, chain_id, position, new_aa):
        self.chain_id = chain_id
        self.position = position
        self.new_aa = new_aa

    def accept_residue(self, residue):
        if residue.id[1] == self.position and residue.parent.id == self.chain_id:
            # Just for demonstration, "mutate" by printing the change
            print(f"Position {self.position}: {residue.resname} → {self.new_aa}")
            return 1
        return 1  # Accept all residues

# Create a copy of the structure for the mutant
mutant_structure = structure.copy()

# In a real application, we would modify the structure to represent the mutation
# For demonstration, we'll just highlight the residue at position 7
mutation_position = 7  # Position in the sequence (1-indexed for PDB)

# Find the residue at position 7 in the first chain
first_chain = next(structure.get_chains())
chain_id = first_chain.id

# Visualization functions
def plot_protein_backbone(structure, title, highlight_pos=None, highlight_color='red'):
    # Extract backbone coordinates (CA atoms)
    ca_atoms = []
    residue_positions = []
    highlighted_ca = None

    for model in structure:
        for chain in model:
            for i, residue in enumerate(chain):
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'].get_coord())
                    residue_positions.append(residue.id[1])
                    if highlight_pos and residue.id[1] == highlight_pos:
                        highlighted_ca = residue['CA'].get_coord()

    # Convert to numpy array for plotting
    if ca_atoms:
        ca_coords = np.array(ca_atoms)

        # Create 3D plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot backbone
        ax.plot(ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2], color='blue', linewidth=2)

        # Highlight specific residue
        if highlighted_ca is not None:
            ax.scatter(highlighted_ca[0], highlighted_ca[1], highlighted_ca[2], 
                      color=highlight_color, s=100, label=f'Position {highlight_pos}')

        ax.set_title(title)
        ax.legend()

        return fig
    else:
        print("No CA atoms found in the structure")
        return None

# Plot wild-type structure
wt_fig = plot_protein_backbone(structure, "Wild-Type GLP-R")
if wt_fig:
    wt_fig.savefig('protein_data/wt_glpr_structure.png')
    plt.close(wt_fig)

# Plot mutant structure with highlighted mutation
mutant_fig = plot_protein_backbone(structure, "P7L Mutant GLP-R", 
                                  highlight_pos=mutation_position, 
                                  highlight_color='red')
if mutant_fig:
    mutant_fig.savefig('protein_data/mutant_glpr_structure.png')
    plt.close(mutant_fig)

# Create a side-by-side comparison
def create_comparison_figure():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7), subplot_kw={'projection': '3d'})

    # Extract CA atoms
    ca_atoms = []
    highlighted_ca = None

    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'].get_coord())
                    if residue.id[1] == mutation_position:
                        highlighted_ca = residue['CA'].get_coord()

    if ca_atoms:
        ca_coords = np.array(ca_atoms)

        # Wild type plot
        ax1.plot(ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2], color='blue', linewidth=2)
        ax1.set_title("Wild-Type GLP-R")

        # Mutant plot
        ax2.plot(ca_coords[:, 0], ca_coords[:, 1], ca_coords[:, 2], color='blue', linewidth=2)
        if highlighted_ca is not None:
            ax2.scatter(highlighted_ca[0], highlighted_ca[1], highlighted_ca[2], 
                       color='red', s=100, label=f'P7L Mutation')
        ax2.set_title("P7L Mutant GLP-R")
        ax2.legend()

        # Set the same view for both plots
        ax1.view_init(elev=30, azim=60)
        ax2.view_init(elev=30, azim=60)

        plt.tight_layout()
        return fig
    else:
        return None

comparison_fig = create_comparison_figure()
if comparison_fig:
    comparison_fig.savefig('protein_data/glpr_comparison.png')
    plt.show()

print("\nVisualization complete!")
print("Three images have been created:")
print("1. protein_data/wt_glpr_structure.png - Wild-type GLP-R structure")
print("2. protein_data/mutant_glpr_structure.png - P7L mutant GLP-R with mutation highlighted")
print("3. protein_data/glpr_comparison.png - Side-by-side comparison")
print("\nIn a real-world protein folding prediction, we would use AlphaFold or similar tools")
print("to predict more accurate structural changes resulting from the mutation.")