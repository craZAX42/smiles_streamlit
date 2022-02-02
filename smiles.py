# pyright: reportMissingImports=false
import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import SessionState
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import hydralit_components as hc
import copy
from rdkit.Chem.Draw import rdDepictor

st.set_page_config(layout='wide',)

menu_data = [
    {'id':'SMILES','icon':"🔬",'label':"SMILES"},
    {'id':'About','icon': "🧑", 'label':"About"},
]

over_theme = {'txc_inactive': '#FFFFFF'}
menu_id = hc.nav_bar(
    menu_definition=menu_data,
    override_theme=over_theme,
    hide_streamlit_markers=False,
    sticky_nav=True,
    sticky_mode='pinned',
)

if menu_id == 'SMILES':
    st.title('Visualize and Compare SMILES in 3D :smile:')
    st.write('View small molecules in 3D straight from a SMILES string. RDKit is used to generate coordinates and py3Dmol for structure rendering. A maximum of 10 conformers for each molecule is generated. To choose different conformers or previously entered molecules, they can be selected from the selectbox on the left. Hydrogen atoms (including polar) are not displayed for clarity.')
    st.write('Only organic molecules are supported. Organometallics will cause an error.')

    ss = SessionState.get(smile_models={}, smile_strings=[], conformers=[])

    HTML_NAME = "viz.html"
    HEIGHT_3D = 500
    WIDTH_3D = 800

    ## 3D MODEL STYLES ##
    styles = ('stick', 'sphere', 'ball and stick')
    pymol_style = st.selectbox('Select 3D Style', styles)

    def generate_conformers(smile):
        """Create a mol and list of conformer IDs from a smile and store it in the Session State.

        Args:
            smile (string): The smile string to generate a molecule and conformers from
        """
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            raise BaseException("Invalid Smile String")
        ss.smile_models[smile] = {'original_mol': mol, 'cids':{}}
        mol = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol)
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
        mol = Chem.RemoveHs(mol)
        ss.smile_models[smile]['conf_mol'] = mol
        for cid in cids:
            ss.smile_models[smile]['cids'][cid] = {}
        
    def create_model(smile, cid, style):
        """Creates a py3Dmol view and writes the result to an HTML file.

        Args:
            m (rdkit molecule): The mol to generate a MolBlock from
            cid (int): the conformer ID to use to generate the MolBlock
            style (string): the style of model to generate
        """
        view = py3Dmol.view(width=WIDTH_3D, height=HEIGHT_3D)
        mb = Chem.MolToMolBlock(ss.smile_models[smile]['conf_mol'], confId=cid)
        view.addModel(mb, 'mol', {'keepH':False})
        if style == 'ball and stick':
            view.setStyle({'sphere':{'radius':0.5}, 'stick':{}})
        elif style == 'sphere':
            view.setStyle({'sphere':{'radius':1}})
        else:
            view.setStyle({style:{}})
        view.zoomTo()
        view.show()
        view.render()
        t = view.js()
        f = open(HTML_NAME, 'w')
        f.write(t.startjs)
        f.write(t.endjs)
        f.close()
        HtmlFile = open(HTML_NAME, 'r', encoding='utf-8')
        source_code = HtmlFile.read()
        ss.smile_models[smile]['cids'][cid][style] = source_code
        
    def add_smile_and_conformers(smile):
        """Take the smile the user input and add it to the list of smile strings
        along with their 2d and 3d structure.

        Args:
            smile (str): The smile string provided by the user
        """
        if smile in ss.smile_strings:
            return
        generate_conformers(smile)
        ss.smile_strings.append(smile)
        
        for cid in ss.smile_models[smile]['cids'].keys():
            for style in styles:
                create_model(smile, cid, style)

    def display_smiles():
        """Generates the 2d and 3d views of the two smiles selected from the sidebar.
        """
        c1,c2=st.columns([1,2])
        with c1:
            AllChem.GenerateDepictionMatching3DStructure(ss.smile_models[first_smile]['original_mol'], ss.smile_models[first_smile]['conf_mol'])
            Draw.MolToFile(ss.smile_models[first_smile]['original_mol'], 'mol.png',ignoreHs=True)
            st.image('mol.png')
        with c2:
            components.html(ss.smile_models[first_smile]['cids'][first_conf_id][pymol_style],width=WIDTH_3D, height=HEIGHT_3D)

        if first_smile == second_smile and first_conf_id == second_conf_id:
            return
        
        ss.smile_models[second_smile]['mol'] = Chem.RemoveHs(ss.smile_models[second_smile]['conf_mol'])
        c1,c2=st.columns([1,2])
        with c1:
            AllChem.GenerateDepictionMatching3DStructure(ss.smile_models[second_smile]['original_mol'], ss.smile_models[second_smile]['conf_mol'])
            Draw.MolToFile(ss.smile_models[second_smile]['original_mol'], 'mol2.png')
            st.image('mol2.png')
        with c2:
            components.html(ss.smile_models[second_smile]['cids'][second_conf_id][pymol_style],width=WIDTH_3D, height=HEIGHT_3D)

    input_smile=st.text_input('Enter SMILES string\nLeft click, hold, then move mouse to rotate 3D view. \
                                Right click, hold then move mouse to zoom in/out or use the scroll wheel.',\
                                'O=C1C2=C(N=CN2C)N(C)C(N1C)=O')

    try:
        add_smile_and_conformers(input_smile)
    except BaseException:
            st.write('Invalid SMILES, please input a valid SMILES string.')
        
    ################## Sidebar ######################
    second_index=[len(ss.smile_strings)-2 if len(ss.smile_strings) > 1 else 0][0]
    first_smile = st.sidebar.selectbox('Select your first desired SMILE', ss.smile_strings, index=len(ss.smile_strings)-1)
    second_smile = st.sidebar.selectbox('Select your second desired SMILE', ss.smile_strings, index=second_index)

    first_conf_id = st.sidebar.selectbox(first_smile + " conformer", ss.smile_models[first_smile]['cids'].keys(), key=1)
    second_conf_id = st.sidebar.selectbox(second_smile + " conformer", ss.smile_models[second_smile]['cids'].keys(), key=2)

    ########## Display 2D and 3D Structures ##########
    display_smiles()

elif menu_id == 'About':
    st.info("CraZAX is a superstar undergrad student in a computational chemistry lab working to build capabilties and automate workflows.")
    st.info("Although CraZAX is majoring in computer science, this project demonstrates his ability to solve an out of domain chemistry problem implementing programmatic solutions. I overlapped as a PhD student in the chemistry lab and now work for a big pharma company. Tech/Pharma companies, get your hands on CraZAX while you can - tenacious-dk.")
