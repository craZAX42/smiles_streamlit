# pyright: reportMissingImports=false
import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import SessionState
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import hydralit_components as hc

st.set_page_config(layout='wide',initial_sidebar_state='collapsed',)

# menu_data = [
#     {'icon': "far fa-copy", 'label':"Left End"},
#     {'id':'Copy','icon':"🐙",'label':"Copy"},
#     {'icon': "fa-solid fa-radar",'label':"Dropdown1", 'submenu':[{'id':' subid11','icon': "fa fa-paperclip", 'label':"Sub-item 1"},{'id':'subid12','icon': "💀", 'label':"Sub-item 2"},{'id':'subid13','icon': "fa fa-database", 'label':"Sub-item 3"}]},
# ]

# over_theme = {'txc_inactive': '#FFFFFF'}
# menu_id = hc.nav_bar(
#     menu_definition=menu_data,
#     override_theme=over_theme,
#     home_name='Home',
#     login_name='Logout',
#     hide_streamlit_markers=False, #will show the st hamburger as well as the navbar now!
#     sticky_nav=True, #at the top or not
#     sticky_mode='pinned', #jumpy or not-jumpy, but sticky or pinned
# )

st.title('SMILES Comparison :smile:')

ss = SessionState.get(smile_models={}, smile_strings=[], conformers=[])

HTML_NAME = "viz.html"

## 3D MODEL STYLES ##
styles = ('stick', 'sphere', 'ball and stick')
pymol_style = st.selectbox('Select 3D Style', styles)

def generate_conformers(smile):
    """Create a mol and list of conformer IDs from a smile and store it in the Session State.

    Args:
        smile (string): The smile string to generate a molecule and conformers from
    """
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol)
    ss.smile_models[smile]['mol'] = mol
    #res = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
    #print("Res:", res)
    for cid in cids:
        ss.smile_models[smile]['cids'][cid] = {}
    
def create_model(m, cid, style):
    """Creates a py3Dmol view and writes the result to an HTML file.

    Args:
        m (rdkit molecule): The mol to generate a MolBlock from
        cid (int): the conformer ID to use to generate the MolBlock
        style (string): the style of model to generate
    """
    view = py3Dmol.view(width=450, height=450)
    mb = Chem.MolToMolBlock(m, confId=cid)
    view.addModel(mb, 'mol', {'keepH':True}) # Set KeepH to False later
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
    
def add_smile_and_conformers(smile):
    """Take the smile the user input and add it to the list of smile strings
    along with their 2d and 3d structure.

    Args:
        smile (str): The smile string provided by the user
    """
    if smile in ss.smile_strings:
        return
    ss.smile_strings.append(smile)
    ss.smile_models[smile] = {'cids':{}}
    generate_conformers(smile)

    for cid in ss.smile_models[smile]['cids'].keys():
        for style in styles:
            create_model(ss.smile_models[smile]['mol'], cid, style)
            HtmlFile = open(HTML_NAME, 'r', encoding='utf-8')
            source_code = HtmlFile.read()
            ss.smile_models[smile]['cids'][cid][style] = source_code
        

def display_smiles():
    """Generates the 2d and 3d views of the two smiles selected from the sidebar.
    """
    c1,c2=st.columns(2)
    with c1:
        Draw.MolToFile(ss.smile_models[first_smile]['mol'], 'mol.png')
        st.image('mol.png')
    with c2:
        components.html(ss.smile_models[first_smile]['cids'][first_conf_id][pymol_style],height=450,width=450)

    if first_smile == second_smile:
        return
    
    c1,c2=st.columns(2)
    with c1:
        Draw.MolToFile(ss.smile_models[second_smile]['mol'], 'mol2.png')
        st.image('mol2.png')
    with c2:
        components.html(ss.smile_models[second_smile]['cids'][second_conf_id][pymol_style], height=450,width=450)

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

first_conf_id = st.sidebar.selectbox(first_smile, ss.smile_models[first_smile]['cids'].keys(), key=1)
second_conf_id = st.sidebar.selectbox(second_smile, ss.smile_models[second_smile]['cids'].keys(), key=2)

# print(ss.smile_strings) # For Debugging
# print(ss.smile_models) # For Debugging

########## Display 2D and 3D Structures ##########
display_smiles()
