
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

adata = sc.read_h5ad('data/BBXB/Tcell_Fibro_'+img+'_spatial.h5ad')

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

# 4. get merge cluster list and selected color from user

st.sidebar.header("Create merged cluster to plot")
st.write('### Plotting merged clusters')

n_clusters = st.sidebar.text_input('How many clusters for new plot',2)

adata_merge_ls = []
color_ls = []
name_ls = []
for i in range(int(n_clusters)):
    cl = st.sidebar.text_input('merged list of cluster '+str(i+1)+ '(using / )', str(i+3)+'/'+str(i+4) )
    name = st.sidebar.text_input('name of merged cluster '+str(i+1), 'Cluster'+str(i+1) )
    color = st.sidebar.color_picker('color merged cluster '+str(i+1), '#ff34ff')
    cl_merge_ls = cl.split('/')
    adata_merge = adata[adata.obs['Annot_level2_label'].isin(cl_merge_ls)]
    adata_merge.obs['merged'] = name
    adata_merge_ls.append(adata_merge)
    color_ls.append(color)
    name_ls.append(name)

adata_mergedall = sc.concat(adata_merge_ls)
adata_mergedall.obs['merged'] = adata_mergedall.obs['merged'].astype('category')
adata_mergedall.obs['merged'] = adata_mergedall.obs['merged'].cat.reorder_categories(name_ls)
adata_mergedall.uns['merged_colors'] = color_ls
fig3 = plt_spatial(adata_mergedall, 'merged')
st.pyplot(fig3)

