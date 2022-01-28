import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import SessionState
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import hydralit_components as hc

#test2
#test
# @st.cache(allow_output_mutation=True)
# def smile_models():
#     return {}

# @st.cache(allow_output_mutation=True, hash_funcs={set: lambda _: None})
# def smile_strings():
#     return []

st.set_page_config(layout='wide',initial_sidebar_state='collapsed',)

menu_data = [
    {'icon': "far fa-copy", 'label':"Left End"},
    {'id':'Copy','icon':"🐙",'label':"Copy"},
    {'icon': "fa-solid fa-radar",'label':"Dropdown1", 'submenu':[{'id':' subid11','icon': "fa fa-paperclip", 'label':"Sub-item 1"},{'id':'subid12','icon': "💀", 'label':"Sub-item 2"},{'id':'subid13','icon': "fa fa-database", 'label':"Sub-item 3"}]},
]

over_theme = {'txc_inactive': '#FFFFFF'}
menu_id = hc.nav_bar(
    menu_definition=menu_data,
    override_theme=over_theme,
    home_name='Home',
    login_name='Logout',
    hide_streamlit_markers=False, #will show the st hamburger as well as the navbar now!
    sticky_nav=True, #at the top or not
    sticky_mode='pinned', #jumpy or not-jumpy, but sticky or pinned
)

st.title('SMILES  + RDKit + Py3DMOL :smiley:')

ss = SessionState.get(smile_models={}, smile_strings=[], conformers=[])

def get_pymol_style(style):
    return style

def generate_HTML_name(style):
    return "viz_" + style + ".html"

styles = ('stick', 'sphere', 'ball and stick')
pymol_style = st.selectbox('Select 3D Style', styles)

def create_models(smi, style):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    cids = Chem.rdDistGeom.EmbedMultipleConfs(mol)
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200, confId=cid)
    ss.smile_models[smi][1]['cids'] = cids
    ###############################################################
    #   =============== TODO ===============
    # * save view in our SessionState
    # * Generate HTML from view when trying to display molecules
    # * Add sidebar options for different conformer selections
    ###############################################################
    mblock = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=450, height=450)
    view.addModel(mblock, 'mol', {'keepH':False})
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
    html_name = generate_HTML_name(style)
    f = open(html_name, 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()
    
def add_smiles(input_smiles):
    """Take the smile the user input and add it to the list of smile strings
    along with their 2d and 3d structure.

    Args:
        input_smiles (str): The string input of a smile provided by the user
    """
    if input_smiles in ss.smile_strings:
        return
    m = Chem.MolFromSmiles(input_smiles)
    ss.smile_strings.append(input_smiles)
    ss.smile_models[input_smiles] = (m, {})
    
    for style in styles:
        create_models(input_smiles, style)
        html_name = generate_HTML_name(style)
        HtmlFile = open(html_name, 'r', encoding='utf-8')
        source_code = HtmlFile.read()
        ss.smile_models[input_smiles][1][style] = source_code
        

def display_smiles(first_smile, second_smile):
    
    c1,c2=st.columns(2)
    with c1:
        Draw.MolToFile(ss.smile_models[first_smile][0], 'mol.png')
        st.image('mol.png')
    with c2:
        components.html(ss.smile_models[first_smile][1][pymol_style],height=450,width=450)
    
    # print("first:", first_smile, "\tsecond:", second_smile)
    
    if first_smile == second_smile:
        return
    
    c1,c2=st.columns(2)
    with c1:
        Draw.MolToFile(ss.smile_models[second_smile][0], 'mol2.png')
        st.image('mol2.png')
    with c2:
        components.html(ss.smile_models[second_smile][1][pymol_style], height=450,width=450)

input_smiles=st.text_input('Enter SMILES string\nLeft click, hold, then move mouse to rotate 3D view. \
                            Right click, hold then move mouse to zoom in/out.',\
                            'O=C1C2=C(N=CN2C)N(C)C(N1C)=O')
# print(input_smiles) # For debugging
try:
    add_smiles(input_smiles)
except BaseException:
        st.write('Invalid SMILES, please input a valid SMILES string.')
    
################ Sidebar ####################
second_index=[len(ss.smile_strings)-2 if len(ss.smile_strings) > 1 else 0][0]
first_smile = st.sidebar.selectbox('Select your first desired SMILE', ss.smile_strings, index=len(ss.smile_strings)-1)
second_smile = st.sidebar.selectbox('Select your second desired SMILE', ss.smile_strings, index=second_index)

first_conf_id = st.sidebar.selectbox(first_smile, ss.smile_models[first_smile][1]['cids'])
if first_smile != second_smile:
    second_conf_id = st.sidebar.selectbox(second_smile, ss.smile_models[second_smile][1]['cids'])

display_smiles(first_smile, second_smile)
