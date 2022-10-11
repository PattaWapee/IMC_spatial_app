from matplotlib import pyplot as plt
import streamlit as st
import anndata as ad
import scanpy as sc
import squidpy as sq
import pandas as pd


st.markdown("""
# Epithelial & Fibroblast subcell classes
"""
            )

st.sidebar.markdown("# Epithelial vs Fibroblast")


def plt_spatial(adata, obs_col):
    '''
    plot spatial graph from adata
    color is labeled from chosen obs column
    '''
    fig = sc.pl.spatial(adata, color=obs_col,
                        title=[obs_col], spot_size=10,
                        show=False, return_fig=True)
    fig = fig[0].get_figure()
    return fig


def filter_adata(adata, obs_col, obs_value_ls):
    '''
    Filtering anndata by list of obs value
    '''
    adata_filtered = adata[adata.obs[obs_col].isin(obs_value_ls)]
    return adata_filtered

# 1. import anndata


img_name = st.sidebar.selectbox(
    "Select image", ('R2(Basal like)', 'R4(ER+)', 'R5(HER2+)'))
img = img_name.split('(')[0]

adata = sc.read_h5ad('data/BBXB/Epi_Fibro_'+img+'_spatial.h5ad')

# 2. choose marker and show spatial plot for all cluster
epi_marker = st.sidebar.selectbox(
    "Select epithelial marker label", adata.obs.columns)
fig = plt_spatial(adata, epi_marker)
st.write('## Spatial plot for ' + img)
st.pyplot(fig)

# 3. get cluster list input from user
cl_input = st.sidebar.text_input('Select clusters for plotting', 'pos,1')
cl_in_ls = cl_input.split(',')

adata_filtered = filter_adata(adata, epi_marker, cl_in_ls)
st.write('### Show only cluster ' + cl_input)
fig2 = plt_spatial(adata_filtered, epi_marker)
st.pyplot(fig2)
# 4. get merge cluster list and selected color from user

st.sidebar.header("Create merged cluster to plot")
st.write('### Plotting merged clusters')

n_clusters = st.sidebar.text_input('How many clusters for new plot', 2)

adata_merge_ls = []
color_ls = []
name_ls = []
for i in range(int(n_clusters)):
    cl = st.sidebar.text_input(
        'merged list of cluster '+str(i+1), str(i+3)+','+str(i+4))
    name = st.sidebar.text_input(
        'name of merged cluster '+str(i+1), 'Cluster'+str(i+1))
    color = st.sidebar.color_picker(
        'color merged cluster '+str(i+1), '#ff34ff')
    cl_merge_ls = cl.split(',')
    adata_merge = filter_adata(adata, epi_marker, cl_merge_ls)
    adata_merge.obs['merged'] = name
    adata_merge_ls.append(adata_merge)
    color_ls.append(color)
    name_ls.append(name)

adata_mergedall = sc.concat(adata_merge_ls)
adata_mergedall.obs['merged'] = adata_mergedall.obs['merged'].astype(
    'category')
adata_mergedall.obs['merged'] = adata_mergedall.obs['merged'].cat.reorder_categories(
    name_ls)
adata_mergedall.uns['merged_colors'] = color_ls
fig3 = plt_spatial(adata_mergedall, 'merged')
st.pyplot(fig3)
