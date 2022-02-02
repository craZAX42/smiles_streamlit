# smiles_streamlit
A streamlit webpage to compare two SMILES strings in 3D.

[Here](https://share.streamlit.io/craZAX42/smiles_streamlit/main/smiles.py) is a link to a Streamlit Webserver running this repo.

If you want to run this on your own local machine, you can clone the repo and run `streamlit run smiles.py` from the directory with the `smiles.py` file.

## How It Works
We're using RDKit to generate 3D coordinates of SMILES strings and py3Dmol renders the structure.

After you input SMILES strings to the textbox, you can compare the 2D and 3D structures side by side. There's a dropdown menu to select the style you want to render the models in. In the sidebar, you can select any two SMILES strings you've already input.

[The input and sidebar] <img src="https://github.com/craZAX42/smiles_streamlit/blob/main/images/navbar_style_input.png"> 

[Two molecules]<img src="https://github.com/craZAX42/smiles_streamlit/blob/main/images/3dview_and_sidebar.png">
