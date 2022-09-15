

from matplotlib import pyplot as plt
import streamlit as st
import anndata as ad
import scanpy as sc
import squidpy as sq

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
