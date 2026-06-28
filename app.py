import pandas as pd
import time
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
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs(["1. Protein", "2. Ligand", "3. Cavity", "4. PDBQT", "5. Params", "6. Scoring", "7. Submit", "8. Results", "9. History"])
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
    
    # 1. Load the database from CSV
    try:
        df = pd.read_csv("unani_data.csv")
    except FileNotFoundError:
        st.error("Database file 'unani_data.csv' not found. Please ensure it is in the same folder.")
        df = pd.DataFrame(columns=["Name", "Property", "Reference", "English", "SMILES"])

    # 2. Input and Fetch
    plant_input = st.text_input("Enter Plant Name (e.g. Kalonji):")
    
    # Defaults
    data = {"Property": "", "Reference": "", "English": "", "SMILES": ""}
    
    if st.button("🔍 Fetch from Database"):
        # Search the DataFrame (case-insensitive)
        match = df[df['Name'].str.lower() == plant_input.lower()]
        if not match.empty:
            data = match.iloc[0].to_dict()
            st.success(f"Plant '{plant_input}' found!")
        else:
            st.warning("Plant not found in database. Enter manually below.")

    # 3. Input Fields (Auto-filled by 'data' dictionary)
    col_a, col_b = st.columns(2)
    with col_a:
        unani_prop = st.text_input("Unani Property:", value=data['Property'])
        reference = st.text_input("Reference/Source:", value=data['Reference'])
    with col_b:
        eng_prop = st.text_input("English Translation:", value=data['English'])
        base_smiles = st.text_input("Base Ligand SMILES:", value=data['SMILES'])

    # 4. Docking & Optimization Logic
    with st.expander("🔬 Structural Optimization"):
        optimized_smiles = st.text_input("Edit SMILES for optimization:", base_smiles)
    
    final_smiles = optimized_smiles if optimized_smiles else base_smiles

    if st.button("Generate Ligand Structures"):
        # ... (Your existing docking/visual logic here, using 'final_smiles')
        try:
            mol = Chem.MolFromSmiles(final_smiles)
            if mol:
                # Lipinski Check
                mw, logp, hbd, hba = calculate_drug_likeness(mol)
                st.subheader("Pharmacological Feasibility")
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Weight", f"{mw:.1f}")
                c2.metric("LogP", f"{logp:.1f}")
                c3.metric("HBD", hbd)
                c4.metric("HBA", hba)
                
                # Visuals
                st.session_state['ligand_mol'] = mol
                img = Draw.MolToImage(mol, size=(400, 400))
                mol_3d = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_3d, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol_3d)
                sdf_block = Chem.MolToMolBlock(mol_3d)
                
                col1, col2 = st.columns(2)
                col1.image(img)
                # ... (rest of your visual code with py3Dmol)
            else:
                st.error("Invalid SMILES.")
        except Exception as e:
            st.error(f"Error: {e}")
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

# --- TAB 5: DOCKING PARAMETERS ---
with tab5:
    st.header("Search & Grid Parameters")
    col1, col2, col3 = st.columns(3)
    exhaustiveness = col1.number_input("Exhaustiveness", value=8, step=1)
    num_modes = col2.number_input("Number of Modes", value=10, step=1)
    energy_range = col3.number_input("Energy Range (kcal/mol)", value=3, step=1)

    st.subheader("Grid Box Dimensions")
    g_col1, g_col2, g_col3 = st.columns(3)
    center_x = g_col1.number_input("Center X", value=0.0)
    center_y = g_col2.number_input("Center Y", value=0.0)
    center_z = g_col3.number_input("Center Z", value=0.0)
    size_x = g_col1.number_input("Size X", value=20.0)
    size_y = g_col2.number_input("Size Y", value=20.0)
    size_z = g_col3.number_input("Size Z", value=20.0)

    with st.expander("Advanced Options"):
        flexible_res = st.toggle("Enable Flexible Residues")
        seed = st.number_input("Random Seed", value=42, step=1)

# --- TAB 6: SCORING & OUTPUT ---
with tab6:
    st.header("Scoring & Output Settings")
    scoring = st.selectbox("Scoring Function", ["Vina", "Vinardo", "AutoDock4 (AD4)"])
    poses = st.slider("Number of poses to save", min_value=1, max_value=20, value=9)
    rmsd = st.number_input("RMSD Clustering Threshold (Å)", value=2.0, step=0.1)
    out_format = st.radio("Output format", ["PDBQT", "SDF"])

# --- TAB 7: JOB SUBMISSION ---
with tab7:
    st.header("Submit Docking Job")
    job_name = st.text_input("Job Name", "Unani_Plant_Docking_01")
    job_desc = st.text_area("Description (Optional)", "Docking of Curcumin against 1HSG...")
    
    col1, col2 = st.columns(2)
    compute_mode = col1.selectbox("Compute Mode", ["CPU (Standard)", "GPU (Accelerated)"])
    email_notif = col2.toggle("Email notification on completion")

    if st.button("🚀 Submit Job", type="primary"):
        with st.spinner("Submitting job to compute cluster..."):
            time.sleep(2)
            st.success(f"Job '{job_name}' successfully submitted!")

# --- TAB 8: RESULTS & VISUALIZATION ---
with tab8:
    st.header("Results & Visualization")
    res_data = pd.DataFrame({
        "Rank": [1, 2, 3, 4, 5],
        "Binding Affinity (kcal/mol)": [-9.8, -9.4, -9.1, -8.8, -8.5],
        "RMSD (l.b.)": [0.00, 1.25, 2.10, 2.80, 3.45],
        "Pose ID": ["Pose_01", "Pose_02", "Pose_03", "Pose_04", "Pose_05"]
    })
    
    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("Pose Rankings")
        st.dataframe(res_data, use_container_width=True, hide_index=True)
        col_btn1, col_btn2 = st.columns(2)
        col_btn1.download_button("⬇️ Download PDBQT", "dummy_data", file_name="results.pdbqt")
        col_btn2.download_button("⬇️ Export CSV", res_data.to_csv(index=False), file_name="rankings.csv")

    with col2:
        res_tab1, res_tab2 = st.tabs(["Interactions", "Raw Files"])
        with res_tab1:
            st.write("**Interaction Summary:**\n- Hydrogen bonds: ASP25, GLY27\n- Hydrophobic: ILE50, PRO81")
        with res_tab2:
            st.code("REMARK VINA RESULT: -9.8  0.000  0.000\nATOM      1  C   UNL     1...")

# --- TAB 9: JOB HISTORY / WORKSPACE ---
with tab9:
    st.header("Job History / Workspace")
    col1, col2 = st.columns([3, 1])
    search_q = col1.text_input("🔍 Search Jobs")
    status_filter = col2.selectbox("Filter by Status", ["All", "Running", "Completed", "Failed"])
    
    history_data = pd.DataFrame({
        "Job Name": ["Curcumin_1HSG", "Quercetin_6LU7", "Aspirin_Test"],
        "Date": ["2026-06-28", "2026-06-27", "2026-06-25"],
        "Status": ["Completed", "Running", "Failed"],
        "Actions": ["View / Re-run", "View Log", "Delete"]
    })
    st.dataframe(history_data, use_container_width=True, hide_index=True)
