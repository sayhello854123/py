#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc
data_type = 'float32'

# this line forces theano to use the GPU and should go before importing cell2location
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=' + data_type + ',force_device=True'

import cell2location

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


# Set paths to data and results used through the document:
sp_data_folder = '/net/data.isilon/ag-saez/bq_shared/scellMI/raw_visium/'
results_folder = '/net/data.isilon/ag-saez/bq_rramirez/visiumMI_revisions/results/deconvolution/c2l/'
regression_model_output = 'rmodelRegressionGeneBackgroundCoverageTorch_26covariates_60574cells_12394genestestMI'
reg_path = f'{results_folder}{regression_model_output}/'
sample_names = ["157771", "157772", "157775", "157777", "157779", "157781", "157782", "157785"]


# In[ ]:


# Reading and analyzing visium data
def read_and_qc(sample_name, path=sp_data_folder + 'rawdata/'):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """

    adata = sc.read_visium(path + str(sample_name),
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']

    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"]                           + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata


# In[ ]:


# Read the data into anndata objects
slides = []
for i in sample_names:
    slides.append(read_and_qc(i, path=sp_data_folder))


# In[ ]:


# Combine anndata objects together
adata = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_names,
    index_unique=None
)

adata.var_names = adata.var['SYMBOL']

adata.var_names_make_unique()

# In[ ]:


# mitochondria-encoded (MT) genes should be removed for spatial mapping
adata.obsm['mt'] = adata[:, adata.var['mt'].values].X.toarray()
adata = adata[:, ~adata.var['mt'].values]


# In[ ]:


adata_vis = adata.copy()
adata_vis.raw = adata_vis


# In[ ]:


## snRNAseq reference (raw counts)
adata_snrna_raw = sc.read(f'{reg_path}sc.h5ad')


# In[ ]:


# Column name containing cell type annotations
covariate_col_names = 'deconv_col'

# Extract a pd.DataFrame with signatures from anndata object
inf_aver = adata_snrna_raw.raw.var.copy()
inf_aver = inf_aver.loc[:, [f'means_cov_effect_{covariate_col_names}_{i}' for i in adata_snrna_raw.obs[covariate_col_names].unique()]]
from re import sub
inf_aver.columns = [sub(f'means_cov_effect_{covariate_col_names}_{i}', '', i) for i in adata_snrna_raw.obs[covariate_col_names].unique()]
inf_aver = inf_aver.iloc[:, inf_aver.columns.argsort()]

# normalise by average experiment scaling factor (corrects for sequencing depth)
inf_aver = inf_aver * adata_snrna_raw.uns['regression_mod']['post_sample_means']['sample_scaling'].mean()


# In[ ]:


# selecting most informative genes based on specificity
selection_specificity = 0.2

# normalise expression signatures:
cell_state_df_norm = (inf_aver.T / inf_aver.sum(1)).T
# apply cut off:
cell_state_df_norm = (cell_state_df_norm > selection_specificity)

# check the number of markers per cell type
cell_state_df_norm.sum(0), (cell_state_df_norm.sum(1) > 0).sum(0)


# In[ ]:


# now we don't need to keep the scRNA-seq data set and a list with slides in memory
del adata_snrna_raw, slides
gc.collect()


# In[ ]:


sc.settings.set_figure_params(dpi = 100, color_map = 'viridis', dpi_save = 100,
                              vector_friendly = True, format = 'pdf',
                              facecolor='white')

r = cell2location.run_cell2location(

      # Single cell reference signatures as pd.DataFrame
      # (could also be data as anndata object for estimating signatures
      #  as cluster average expression - `sc_data=adata_snrna_raw`)
      sc_data=inf_aver,
      # Spatial data as anndata object
      sp_data=adata_vis,

      # the column in sc_data.obs that gives cluster idenitity of each cell
      summ_sc_data_args={'cluster_col': "deconv_col",
                         # select marker genes of cell types by specificity of their expression signatures
                         'selection': "cluster_specificity",
                         # specificity cutoff (1 = max, 0 = min)
                         'selection_specificity': 0.2
                        },

      train_args={'use_raw': True, # By default uses raw slots in both of the input datasets.
                  'n_iter': 40000, # Increase the number of iterations if needed (see QC below)

                  # When analysing the data that contains multiple experiments,
                  # cell2location automatically enters the mode which pools information across experiments
                  'sample_name_col': 'sample'}, # Column in sp_data.obs with experiment ID (see above)


      export_args={'path': results_folder, # path where to save results
                   'run_name_suffix': '' # optinal suffix to modify the name the run
                  },

      model_kwargs={ # Prior on the number of cells, cell types and co-located groups

                    'cell_number_prior': {
                        # - N - the expected number of cells per location:
                        'cells_per_spot': 8,
                        # - A - the expected number of cell types per location:
                        'factors_per_spot': 4,
                        # - Y - the expected number of co-located cell type groups per location
                        'combs_per_spot': 3
                    },

                     # Prior beliefs on the sensitivity of spatial technology:
                    'gene_level_prior':{
                        # Prior on the mean
                        'mean': 1/2,
                        # Prior on standard deviation,
                        # a good choice of this value should be at least 2 times lower that the mean
                        'sd': 1/4
                    }
      }
)

