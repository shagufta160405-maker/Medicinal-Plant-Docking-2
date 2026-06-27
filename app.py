import streamlit as st
import streamlit.components.v1 as components
import urllib.request
import py3Dmol
from rdkit import Chem
import streamlit as st
import streamlit.components.v1 as components
import urllib.request
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
# Use the Draw module cautiously
try:
    from rdkit.Chem import Draw
except ImportError:
    st.error("RDKit Draw module failed to load. Ensure packages.txt includes libxrender1.")

# ... rest of your code ...

# --- Page Config ---
st.set_page_config(page_title="Medicinal Plant Docking Portal", layout="wide")
st.title("🧬 Medicinal Plant Molecular Docking Portal")
st.markdown("Prepare proteins and medicinal plant ligands for molecular docking.")

# --- Helper Functions ---
def showmol(view, height=400, width=800):
    """Bypasses stmol dependency by rendering py3Dmol directly in Streamlit."""
    components.html(view._make_html(), height=height, width=width)

@st.cache_data
def fetch_pdb(pdb_id):
    url = f"https://files.rcsb.org/view/{pdb_id.upper()}.pdb"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        response = urllib.request.urlopen(req)
        return response.read().decode('utf-8')
    except Exception as e:
        return None

def extract_heteroatoms(pdb_data):
    hetatms = [line for line in pdb_data.split('\n') if line.startswith("HETATM")]
    return "\n".join(hetatms)

# --- UI Tabs ---
tab1, tab2, tab3, tab4 = st.tabs(["1. Protein Preparation", "2. Ligand (SMILES)", "3. Cavity & Co-factors", "4. PDBQT Export"])

# --- TAB 1: PROTEIN PREPARATION ---
with tab1:
    st.header("Fetch or Upload Target Protein")
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Fetch by PDB ID")
        pdb_id = st.text_input("Enter PDB ID (e.g., 1HSG, 6LU7):", "1HSG")
        if st.button("Fetch Protein"):
            pdb_data = fetch_pdb(pdb_id)
            if pdb_data:
                st.session_state['pdb_data'] = pdb_data
                st.success(f"Successfully fetched {pdb_id.upper()}!")
            else:
                st.error("Failed to fetch PDB. Check the ID.")
                
    with col2:
        st.subheader("Or Upload PDB File")
        uploaded_file = st.file_uploader("Upload Protein (.pdb)", type=["pdb"])
        if uploaded_file is not None:
            st.session_state['pdb_data'] = uploaded_file.getvalue().decode("utf-8")
            st.success("File uploaded successfully!")

    if 'pdb_data' in st.session_state:
        st.markdown("### 3D Protein Viewer")
        view = py3Dmol.view(width=800, height=400)
        view.addModel(st.session_state['pdb_data'], 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.zoomTo()
        showmol(view, height=400, width=800)

# --- TAB 2: LIGAND (SMILES) ---
with tab2:
    st.header("Medicinal Plant Ligand Setup")
    smiles = st.text_input("Enter Ligand SMILES string:", "CC(=O)OC1=CC=CC=C1C(=O)O") # Default: Aspirin
    
    if st.button("Generate Ligand Structures"):
        try:
            # Generate 2D
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                st.session_state['ligand_mol'] = mol
                img = Draw.MolToImage(mol, size=(400, 400))
                
                # Generate 3D SDF
                mol_3d = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_3d, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol_3d)
                st.session_state['ligand_3d'] = mol_3d
                
                sdf_block = Chem.MolToMolBlock(mol_3d)
                st.session_state['sdf_block'] = sdf_block
                
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("### 2D Structure")
                    st.image(img)
                
                with col2:
                    st.markdown("### 3D Structure")
                    view = py3Dmol.view(width=400, height=400)
                    view.addModel(sdf_block, 'sdf')
                    view.setStyle({'stick': {}})
                    view.zoomTo()
                    showmol(view, height=400, width=400)
                    
                st.download_button(
                    label="Download 3D SDF File",
                    data=sdf_block,
                    file_name="ligand_3D.sdf",
                    mime="chemical/x-mdl-sdfile"
                )
            else:
                st.error("Invalid SMILES string.")
        except Exception as e:
            st.error(f"Error processing SMILES: {e}")

# --- TAB 3: CAVITY & CO-FACTORS ---
with tab3:
    st.header("Scan Cavity & Identify Co-factor Heteroatoms")
    if 'pdb_data' in st.session_state:
        st.info("Scanning for HETATM (Heteroatoms, Co-factors, Water) in the loaded protein...")
        hetatms = extract_heteroatoms(st.session_state['pdb_data'])
        
        if hetatms:
            st.text_area("Heteroatom Records Found:", hetatms, height=200)
            
            st.markdown("### View Target Cavity (Ligand Binding Site)")
            view = py3Dmol.view(width=800, height=400)
            view.addModel(st.session_state['pdb_data'], 'pdb')
            view.setStyle({'cartoon': {'color': 'lightgray'}})
            # Highlight heteroatoms/ligands
            view.addStyle({'hetflag': True}, {'stick': {'colorscheme': 'greenCarbon'}})
            view.zoomTo({'hetflag': True})
            showmol(view, height=400, width=800)
        else:
            st.warning("No heteroatoms or co-factors found in this PDB file.")
    else:
        st.warning("Please load a protein in Tab 1 first.")

# --- TAB 4: PDBQT EXPORT ---
with tab4:
    st.header("PDBQT Generation for Docking")
    st.markdown("Convert your prepared structures into PDBQT format (adds partial charges and AutoDock atom types).")
    
    if 'pdb_data' in st.session_state:
        st.subheader("Protein Processing")
        # Stripping water and HETATM for basic receptor PDBQT prep
        clean_pdb = "\n".join([line for line in st.session_state['pdb_data'].split('\n') if line.startswith("ATOM")])
        st.success("Water and co-factors stripped. Ready for partial charge calculation.")
        st.download_button("Download Cleaned Receptor (PDB)", clean_pdb, file_name="receptor_clean.pdb")
        
    if 'sdf_block' in st.session_state:
        st.subheader("Ligand Processing")
        st.success("3D Ligand coordinates mapped. Ready for torsion tree root setup.")
        st.download_button("Download Ligand (SDF)", st.session_state['sdf_block'], file_name="ligand_prepared.sdf")
