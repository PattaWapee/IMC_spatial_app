from matplotlib import pyplot as plt
import streamlit as st
import anndata as ad
import scanpy as sc
import squidpy as sq


st.markdown("""
# Epithelial & Fibroblast subcell classes
"""
            )

sd.sidebar.markdown("# Epithelial vs Fibroblast")


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

adata = sc.read_h5ad('data/Epi_Fibro_'+img+'_spatial.h5ad')

# 2. choose marker and show spatial plot for all cluster
epi_marker = st.sidebar.selectbox(
    "Select epithelial marker label", adata.obs.columns)
fig = plt_spatial(adata, epi_marker)
st.write('## Spatial plot for ' + img)
st.pyplot(fig)

# 3. get cluster list input from user
cl_input = st.sidebar.text_input('Choose cluster for plotting', 'pos,1')
cl_in_ls = cl_input.split(',')

adata_filtered = filter_adata(adata, epi_marker, cl_in_ls)
st.write('### Show only cluster ' + cl_input)
fig2 = plt_spatial(adata_filtered, epi_marker)
st.pyplot(fig2)

