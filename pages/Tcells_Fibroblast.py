
from matplotlib import pyplot as plt
import streamlit as st
import anndata as ad
import scanpy as sc
import squidpy as sq

import sys
sys.path.append('.')
from pages.utils.utils import *

st.markdown("""
# Tcells & Fibroblast subcell classes
"""
            )

# 1. import anndata

img_name = st.sidebar.selectbox(
    "Select image", ('R2(Basal like)', 'R4(ER+)', 'R5(HER2+)'))
img = img_name.split('(')[0]

adata = sc.read_h5ad('data/Tcell_Fibro_'+img+'_spatial.h5ad')

# 2. choose marker and show spatial plot for all cluster
fig = plt_spatial(adata, 'Annot_level2_label')
st.write('## Spatial plot for ' + img)
st.pyplot(fig)


# 3. get cluster list input from user

cl_input = st.sidebar.text_input(
    'Choose cluster for plotting (using / for calling another cluster)', 'CD4+, CD8+/1')
cl_in_ls = cl_input.split('/')

adata_filtered = adata[adata.obs['Annot_level2_label'].isin(cl_in_ls)]
st.write('### Show only cluster ' + cl_input)
fig2 = plt_spatial(adata_filtered, 'Annot_level2_label')
st.pyplot(fig2)
