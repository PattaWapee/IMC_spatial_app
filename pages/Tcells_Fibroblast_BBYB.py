from matplotlib import pyplot as plt
import streamlit as st
import anndata as ad
import scanpy as sc
import squidpy as sq

import sys
sys.path.append('.')
from pages.utils.utils import *

st.markdown("""
# Tcells & Fibroblast BBYB
"""
            )
# 1. import anndata
st.sidebar.header("Section 1. Spatial plot")

img_name = st.sidebar.selectbox(
    "1. Select image", ('R2(Basal like)', 'R4(ER+)', 'R5(HER2+)'))
img = img_name.split('(')[0]

# 2. choose Fibroblast marker type
Fibro_marker = st.sidebar.selectbox(
    "Select Fibroblast marker label", ('aSMA','FAP','PDGFRb'))

adata = sc.read_h5ad('data/BBYB/Tcell_Fibro_'+
                     Fibro_marker + '_'+
                     img+'_spatial.h5ad')

# 2. Show spatial plot for all cluster
fig = plt_spatial(adata, 'Annot_level2_label')
st.write('## 1. Spatial plot for ' + img+ 
         ' Tcells vs Fibroblast ' + Fibro_marker)
st.pyplot(fig)

# 3. get cluster list input from user
st.sidebar.header("Section 2. Choose cluster to plot")

cl_input = st.sidebar.text_input(
    '2. Choose cluster for plotting (using / for calling another cluster)', 'CD4+, CD8+/1')
cl_in_ls = cl_input.split('/')

adata_filtered = adata[adata.obs['Annot_level2_label'].isin(cl_in_ls)]
st.write('## 2. Show only cluster ' + cl_input)
fig2 = plt_spatial(adata_filtered, 'Annot_level2_label')
st.pyplot(fig2)

# 4. get merge cluster list and selected color from user
st.sidebar.header("Section 3. Create merged cluster with manual color")

st.write('## 3. Plotting merged manual clusters')

n_clusters = st.sidebar.text_input('How many clusters for new plot',2)
colors = ['#ffff00', '#1ce6ff', '#ff34ff', '#ff4a46', '#008941', '#006fa6',
       '#a30059', '#ffdbe5', '#7a4900', '#0000a6', '#63ffac', '#b79762',
       '#004d43', '#8fb0ff', '#997d87', '#5a0007', '#809693', '#6a3a4c',
       '#1b4400', '#4fc601', '#3b5dff', '#4a3b53', '#ff2f80', '#61615a',
       '#ba0900', '#6b7900', '#00c2a0', '#ffaa92', '#ff90c9', '#b903aa',
       '#0000FF', '#ff0000']

adata_merge_ls = []
color_ls = []
name_ls = []
try:
    for i in range(int(n_clusters)):
        cl = st.sidebar.text_input('merged list of cluster '+str(i+1)+ '(using / )', str(i+3)+'/'+str(i+4) )
        name = st.sidebar.text_input('name of merged cluster '+str(i+1), 'Cluster'+str(i+1) )
        color = st.sidebar.color_picker('color merged cluster '+str(i+1), colors[i])
        cl_merge_ls = cl.split('/')
        adata_merge = adata[adata.obs['Annot_level2_label'].isin(cl_merge_ls)]
        adata_merge.obs['merged'] = name
        adata_merge_ls.append(adata_merge)
        color_ls.append(color)
        name_ls.append(name)
except:
    st.write('Please select number of cluster less than 33 ')

adata_mergedall = sc.concat(adata_merge_ls)
adata_mergedall.obs['merged'] = adata_mergedall.obs['merged'].astype('category')
adata_mergedall.obs['merged'] = adata_mergedall.obs['merged'].cat.reorder_categories(name_ls)
adata_mergedall.uns['merged_colors'] = color_ls
fig3 = plt_spatial(adata_mergedall, 'merged')
st.pyplot(fig3)

    