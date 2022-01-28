# pyright: reportMissingImports=false
from pages import pages

import streamlit as st
from streamlit_multipage import MultiPage

app = MultiPage()
app.st = st

for app_name, app_tuple in pages.items():
    app.add_app(app_name, app_tuple[0], initial_page=app_tuple[1])
    
app.run()