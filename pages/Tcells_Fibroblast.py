
from matplotlib import pyplot as plt
import streamlit as st
import anndata as ad
import scanpy as sc
import squidpy as sq

from utils.utils import *

st.markdown("""
# Tcells & Fibroblast subcell classes
"""
            )

# 1. import anndata

img_name = st.sidebar.selectbox(
    "Select image", ('R2(Basal like)', 'R4(ER+)', 'R5(HER2+)'))
img = img_name.split('(')[0]

adata = sc.read_h5ad('Tcells/Epi_Fibro_'+img+'_spatial.h5ad')
