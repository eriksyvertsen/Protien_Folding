import streamlit as st
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBList import PDBList
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import PPBuilder, protein_letters_3to1
from Bio.Data.IUPACData import protein_letters_1to3

# Set page configuration
st.set_page_config(
    page_title="GLPR Variant Visualizer",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Create directory for files
if not os.path.exists('protein_data'):
    os.makedirs('protein_data')

# Add title and description
st.title("GLP-R Variant Visualizer")
st.markdown("""
This tool allows you to visualize the location of missense mutations in the GLP-R protein structure.
Enter the variant details below to see where the mutation is located in the 3D structure.
""")

# Function to download a PDB structure 
# We'll use GLP1R structure as example (PDB ID: 6X1A)
@st.cache_resource
def get_structure(pdb_id, file_path):
    pdb_file = os.path.join(file_path, f"pdb{pdb_id.lower()}.ent")

    # Check if the file already exists
    if not os.path.exists(pdb_file):
        st.info(f"Downloading structure for PDB ID: {pdb_id}...")
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(pdb_id, pdir=file_path, file_format='pdb')

    return pdb_file

# Parse the PDB file and return structure
@st.cache_resource
def parse_structure(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('GLPR', pdb_file)
    return structure

# Extract sequence from structure
@st.cache_data
def get_sequence(structure):
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        sequence = pp.get_sequence()
        break  # Just take the first chain's sequence
    return sequence

# Create a mapping of residue positions to structure positions
@st.cache_data
def map_residue_positions(structure):
    residue_map = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':  # Standard amino acid
                    residue_map[residue.get_id()[1]] = {
                        'chain': chain.get_id(),
                        'resname': residue.get_resname(),
                        'coords': residue['CA'].get_coord() if 'CA' in residue else None
                    }
    return residue_map

# Visualization function
def plot_protein_with_mutation(structure, mutation_pos, original_aa, new_aa, title):
    # Extract backbone coordinates (CA atoms)
    ca_atoms = []
    highlighted_ca = None

    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue and residue.get_id()[0] == ' ':  # Standard amino acid
                    ca_atoms.append(residue['CA'].get_coord())
                    if residue.get_id()[1] == mutation_pos:
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
                      color='red', s=100, 
                      label=f'Position {mutation_pos}: {original_aa} → {new_aa}')

        ax.set_title(title)
        ax.legend()

        # Set view angle
        ax.view_init(elev=30, azim=60)

        return fig
    else:
        st.error("No CA atoms found in the structure")
        return None

# Main function to generate visualization
def generate_visualization(position, new_aa):
    # Download structure
    pdb_file = get_structure('6X1A', 'protein_data')

    # Parse structure
    structure = parse_structure(pdb_file)

    # Get sequence
    sequence = get_sequence(structure)

    # Map residue positions
    residue_map = map_residue_positions(structure)

    # Check if the position is valid
    if position not in residue_map:
        st.error(f"Position {position} is not found in the structure. Valid positions are between {min(residue_map.keys())} and {max(residue_map.keys())}.")
        return

    # Get original amino acid
    try:
        original_aa_3letter = residue_map[position]['resname']
        original_aa = protein_letters_3to1[original_aa_3letter]
    except:
        st.error(f"Could not determine the original amino acid at position {position}.")
        return

    # Create visualization
    fig = plot_protein_with_mutation(structure, position, original_aa, new_aa, 
                                    f"GLP-R with {original_aa}{position}{new_aa} Mutation")

    return fig, original_aa, sequence

# Sidebar for inputs
st.sidebar.header("Mutation Parameters")

# Position input
position = st.sidebar.number_input("Position in the sequence", min_value=1, value=7)

# Amino acid selection
amino_acids = list(protein_letters_1to3.keys())
new_aa = st.sidebar.selectbox("New amino acid", amino_acids)

# Generate button
if st.sidebar.button("Generate Visualization"):
    with st.spinner("Generating visualization..."):
        result = generate_visualization(position, new_aa)

        if result:
            fig, original_aa, sequence = result

            # Display mutation info
            st.subheader("Mutation Details")
            st.markdown(f"**Original amino acid:** {original_aa}")
            st.markdown(f"**New amino acid:** {new_aa}")
            st.markdown(f"**Position:** {position}")
            st.markdown(f"**Mutation notation:** {original_aa}{position}{new_aa}")

            # Display figure
            st.pyplot(fig)

            # Display sequence context
            st.subheader("Sequence Context")

            # Determine the context range (10 residues before and after)
            start = max(0, position - 11)
            end = min(len(sequence), position + 10) 

            context_seq = str(sequence)[start:end]
            highlight_pos = position - start - 1

            # Format the sequence with the mutated position highlighted
            formatted_seq = ""
            for i, aa in enumerate(context_seq):
                if i == highlight_pos:
                    formatted_seq += f"<span style='color:red; font-weight:bold;'>{aa}</span>"
                else:
                    formatted_seq += aa

            st.markdown(f"Position {start+1}-{end}:")
            st.markdown(f"<pre>{formatted_seq}</pre>", unsafe_allow_html=True)

# Add information about the structure
st.sidebar.markdown("---")
st.sidebar.subheader("About the Structure")
st.sidebar.markdown("""
This visualization uses the GLP1R structure (PDB ID: 6X1A).

The Glucagon-like peptide 1 receptor (GLP1R) is an important target for the treatment of type 2 diabetes.

Reference: Zhao, P., et al. (2020). Activation of the GLP-1 receptor by a non-peptidic agonist. Nature, 577(7790), 432-436.
""")

# Add GitHub link
st.sidebar.markdown("---")
st.sidebar.markdown("Created with Streamlit, BioPython, and Matplotlib")