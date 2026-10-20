### This script contains the entire CellOracle workflow, adapted to process one cell type based on the SLURM array index

# Python
# Don't change anything here
import os
import sys
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

# Set the working dir
# Change to your own directory
# NOTE: lab cluster path below -- change it to your own path
os.chdir("/resnick/groups/NolanLab/cchen/Feronia/analysis/")

# Import and check CellOracle
# Don't change anything here
import celloracle as co
co.__version__

# Don't change anything here
# Read in an input needed by CellOracle
# NOTE: lab cluster path below -- change it to your own path
id2name = pd.read_table('/resnick/groups/NolanLab/shared_resources/source_files/celloracle_inputs/features.tsv.gz', 
                        compression='gzip', header=None)

# Read in the h5ad object
# Change it to your object name
ADATA = sc.read_h5ad('CellOracle/CellOracle_input/20250509_4S_feronia_salinity_time_course_untreated_celloracle_input.h5ad')
# Load TF info
# Don't change anything here
# NOTE: lab cluster path below -- change it to your own path
wanted_TFs = pd.read_csv("/resnick/groups/NolanLab/shared_resources/source_files/celloracle_inputs/Kay_TF_thalemine_annotations.csv")

# Load base_GRN
# Don't change anything here
import pickle

# NOTE: lab cluster path below -- change it to your own path
with open("/resnick/groups/NolanLab/shared_resources/source_files/celloracle_inputs/base_GRN_TG_to_TF_dictionary.pkl", 'rb') as file_handle:
    TG_to_TF_dictionary = pickle.load(file_handle)
    
DE = pd.read_csv("./muscat_result/20260317_salinity_timecourse_useDateT_noSuspiciousClusters/group_timecelltypeannoLi/fer_4-WT/sig_de_fil.csv", 
                 index_col=None)


# cortex
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
cortex_names = vc[vc.index.str.contains("cortex", case=False, na=False)]
cortex_names.index.tolist()

# endodermis
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
endodermis_names = vc[vc.index.str.contains("endodermis", case=False, na=False)]
endodermis_names.index.tolist()

# trichoblast (need to exclude atrichoblast)
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
trichoblast_names = vc[vc.index.str.contains("trichoblast", case=False, na=False) & ~vc.index.str.contains("atrichoblast", case=False, na=False)]
trichoblast_names.index.tolist()

# atrichoblast
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
atrichoblast_names = vc[vc.index.str.contains("atrichoblast", case=False, na=False)]
atrichoblast_names.index.tolist()

# pericycle
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
pericycle_names = vc[vc.index.str.contains("pericycle", case=False, na=False)]
pericycle_names.index.tolist()

#phloem  (need to exclude pericycle)
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
phloem_names = vc[vc.index.str.contains("phloem", case=False, na=False) & ~vc.index.str.contains("pericycle", case=False, na=False)]
phloem_names.index.tolist()

#xylem (need to exclude pericycle)
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
xylem_names = vc[vc.index.str.contains("xylem", case=False, na=False) & ~vc.index.str.contains("pericycle", case=False, na=False)]
xylem_names.index.tolist()

# procambium
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
procambium_names = vc[vc.index.str.contains("procambium", case=False, na=False)]
procambium_names.index.tolist()

# lrc
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
lrc_names = vc[vc.index.str.contains("Lateral Root Cap", case=False, na=False)]
lrc_names.index.tolist()

# columella
# Don't change anything here
vc = pd.Series(ADATA.obs['time.celltype.anno.Li']).value_counts()
columella_names = vc[vc.index.str.contains("columella", case=False, na=False)]
columella_names.index.tolist()

# List the celltype.time.anno with only 1 count
# They have to be removed manualy in the cell right below because they'll cause error in CellOracle
# Don't change anything here
vc = ADATA.obs['celltype.time.anno'].value_counts()
vc[vc == 1]

# Important: You need to manually exclude (if any) the celltype.time.anno with only 1 count shown above!
# Examples (with ~) are in pericycle and procambium
# Change the lines with "~" (which is the NOT operator) to fit your samples. Remove them if there's nothing to exclude.
# 10 cell types in total

cortex = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(cortex_names.index.tolist()))]

endodermis = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(endodermis_names.index.tolist()))]

trichoblast = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(trichoblast_names.index.tolist()))]

atrichoblast = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(atrichoblast_names.index.tolist()))]

pericycle = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(pericycle_names.index.tolist()))]

phloem = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(phloem_names.index.tolist()))]

xylem = ADATA.obs.index[
    (~ADATA.obs['celltype.time.anno'].isin(['Proliferation Domain_Protoxylem_WT_untreated'])) & 
    (ADATA.obs['time.celltype.anno.Li'].isin(xylem_names.index.tolist()))]

procambium = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(procambium_names.index.tolist()))]

lrc = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(lrc_names.index.tolist()))]

columella = ADATA.obs.index[
    (ADATA.obs['time.celltype.anno.Li'].isin(columella_names.index.tolist()))]

# Build a dictionary to store the cells (barcodes) in each cell type
# Important: include only the cell types you want to analyze! Comment out or remove the unwanted cell types
celltypes = {
    'cortex': cortex,
    #'endodermis': endodermis,
    #'trichoblast': trichoblast,
    'atrichoblast': atrichoblast,
    #'pericycle': pericycle, 
    #'phloem': phloem,
    #'xylem': xylem, 
    #'procambium': procambium, 
    #'lrc': lrc, 
    'columella': columella
}

# Define the directory where all the output will be stored
# NOTE: lab cluster path below -- change it to your own path
output_dir = '/resnick/groups/NolanLab/cchen/Feronia/analysis/CellOracle/CellOracle_output/20260317_feronia_salinity_timecourse_WT_untreated_fer4_untreated_newIntegratedSeurat3_NoSuspiciousCluster'
os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn’t exist
os.chdir(output_dir)  # Change the current working directory

# Get cell type from SLURM array index
cell_type_idx = int(sys.argv[1])
x = list(celltypes.keys())[cell_type_idx]

# Runs CellOracle for the specified cell type
# Don't change anything here
print(f"Processing {x}")

# First, subset the DE enes based on the current cluster/cell type
# Special cases for lrc and trichoblast, so have to process them separately 
if x == 'lrc':
    # Search for 'lateral root cap' in the cluster
    my_cluster = DE['cluster_id'][DE['cluster_id'].str.contains("lateral root cap", case=False, na=False)].unique()
    
elif x == 'trichoblast':
    # Exclude any clusters containing 'atrichoblast' for 'trichoblast'
    my_cluster = DE['cluster_id'][DE['cluster_id'].str.contains(x, case=False, na=False)].unique()
    my_cluster = [cluster for cluster in my_cluster if 'atrichoblast' not in cluster.lower()]
    
else:
    # All other celltypes
    my_cluster = DE['cluster_id'][DE['cluster_id'].str.contains(x, case=False, na=False)].unique()

print(f"Clusters for {x}: {my_cluster}")

# Subset the genes based on the cluster_id (cell type)
gene_of_interest = ADATA.var.index[ADATA.var.index.isin(wanted_TFs['GeneID']) |
                                  ADATA.var.index.isin(DE[DE['cluster_id'].isin(my_cluster)]['gene'].unique())]

# Subset the genes
# Note that here I don't subset genes based on the cluster_id (cell type)! All cell types will use the same set of DE genes!
#gene_of_interest = ADATA.var.index[ADATA.var.index.isin(wanted_TFs['GeneID']) |  
#                               ADATA.var.index.isin(DE['gene'].unique())]  

print('Number of genes : '+format(len(gene_of_interest)))

# Subset the original AnnData object (ADATA) to include only cells of type 'x' and selected genes  
adata = ADATA[celltypes[x], gene_of_interest]
print(adata)

# Perform Principal Component Analysis (PCA) to reduce dimensionality
sc.pp.pca(adata)
# Compute the k-nearest neighbors (KNN) graph for downstream clustering and visualization  
sc.pp.neighbors(adata)
# Compute the UMAP
sc.tl.umap(adata)
# Set the figure size for the UMAP plot  
plt.rcParams["figure.figsize"] = [6, 6]
# Plot the UMAP, coloring cells based on 'celltype.time.anno', with the legend displayed on the plot  
sc.pl.umap(adata, color=['celltype.time.anno'], legend_loc='on data', legend_fontsize='small')
# Rename the '_index' column in the raw AnnData object to 'features'  
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
# Save the processed AnnData object to an HDF5 file with a name based on the cell type 'x'  
adata.write_h5ad(x+'_rawRNA_celloracle.h5ad')

print('Instantiate Oracle object')
# Create an instance of the Oracle class from CellOracle
oracle = co.Oracle()

# In this notebook, we use raw RNA values as an input of Oracle object
# Assign the raw RNA count values to `adata.X` from the 'raw_count' layer  
# This ensures that Oracle receives the unnormalized count data for analysis  
adata.X = adata.layers["raw_count"].copy()

print('Instantiate Oracle object')
# Import the AnnData object into the Oracle instance  
# - `cluster_column_name="celltype.time.anno"` specifies the metadata column containing cell type annotations  
# - `embedding_name="X_umap"` tells Oracle to use the UMAP embedding for visualization  
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="celltype.time.anno",embedding_name="X_umap")

# Import the TF regulatory network data into Oracle  
# This provides a dictionary mapping target genes (TG) to transcription factors (TF)  
oracle.import_TF_data(TFdict=TG_to_TF_dictionary)

print('KNN computation')
# Perform Principal Component Analysis (PCA) on the Oracle object 
oracle.perform_PCA()
# Determine the optimal number of principal components to use  
# - `np.cumsum(oracle.pca.explained_variance_ratio_)` computes the cumulative explained variance  
# - `np.diff()` calculates the change in variance between consecutive components  
# - The condition `> 0.002` finds the point where the increase in explained variance becomes insignificant  
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
# Set an upper limit of 50 principal components to prevent overfitting
n_comps = min(n_comps, 50)
# Get the total number of cells in the dataset
n_cell = oracle.adata.shape[0]
# Print the total number of cells  
print(f"cell number is :{n_cell}")
# Determine the value of `k` for KNN based on 2.5% of the total cell count 
k = int(0.025*n_cell)
# Print the automatically selected k value  
print(f"Auto-selected k is :{k}")
# Perform KNN imputation to reconstruct missing values and refine the data structure  
# - `n_pca_dims=n_comps`: Use the determined number of PCA components  
# - `k=k`: Set the number of neighbors based on the calculated value  
# - `balanced=True`: Ensures balanced neighborhood selection  
# - `b_sight=k*8`: Defines the search space for nearest neighbors  
# - `b_maxl=k*4`: Limits the number of connections  
# - `n_jobs=8`: Utilizes four CPU cores for parallel computation  
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=20)

print('GRN construction')
# Construct gene regulatory network (GRN) using the Oracle object
# Generate regulatory links for GRN based on cell-type annotations.  
# 'alpha' controls the sparsity of the network, 'verbose_level' determines the logging detail, and 'test_mode' runs a test if True.  
links = oracle.get_links(cluster_name_for_GRN_unit="celltype.time.anno", alpha=10, verbose_level=10, test_mode=False)

print('Save Links object')
# Save the links object as an HDF5 file for later use.  
links.to_hdf5(file_path=x+"_rawRNA.links.celloracle.links")

print('Filter Links')
# Filter links based on statistical significance (p-value < 0.01),  
# absolute coefficient weight, and keep only the top 20,000 strongest links.  
links.filter_links(p=0.01, weight="coef_abs", threshold_number=20000)

print('Calculate network scores. It takes several minutes')
# Calculate scores for each regulatory interaction in the network.
links.get_score()
# Display the first few rows of the merged score table to check results.
links.merged_score.head()

print('Save Links object')
# Save the filtered links with a threshold applied  
links.to_hdf5(file_path=x+"_rawRNA_p01_20k.links.celloracle.links")

print('Save merged score')
# Extract and save the final gene regulatory network scores. Store the merged scores in res
res = links.merged_score

# Initialize an empty list to store gene names.  
name = []
# Convert gene IDs into gene names using the id2name mapping.  
for i in res.index.values:
    name.append(id2name[1][np.where(id2name[0]==i)[0]].tolist()[0])
# Insert the gene names as the first column in the results dataframe.  
res.insert(0, 'name', name)
# Save the final gene regulatory network scores as a CSV file.  
res.to_csv(x+"_rawRNA_celloracle_gene_score_iGRN.csv")

# Initialize an empty DataFrame to store network data. 
df = pd.DataFrame()
# Initialize an empty list to store cluster labels.
cluster = []

print('Save network')
for i in links.links_dict.keys():
    # Append the filtered links for each cluster to the DataFrame. 
    df = df.append(links.filtered_links[i])
    # Extend the cluster list with repeated values of the current cluster name,  
    # ensuring each row gets assigned the correct cluster label (assuming 20,000 links per cluster).  
    cluster.extend(np.repeat(i, 20000).tolist())

# Assign the cluster labels to the DataFrame.
df['cluster'] = cluster
# Reorder and select relevant columns for the final output.  
df = df[['source', 'target', 'cluster', 'coef_mean', 'coef_abs', 'p', '-logp']]

# Initialize an empty list to store source gene names.  
name = []
# Convert source gene IDs into gene names using the id2name mapping
for i in df['source']:
    name.append(id2name[1][np.where(id2name[0]==i)[0]].tolist()[0])
# Insert the source gene names as the second column in the DataFrame. 
df.insert(1, 'source_name', name)

# Reinitialize the list to store target gene names.  
name = []
# Convert target gene IDs into gene names using the id2name mapping.
for i in df['target']:
    name.append(id2name[1][np.where(id2name[0]==i)[0]].tolist()[0])
# Insert the target gene names as the fourth column in the DataFrame. 
df.insert(3, 'target_name', name)
# Save the final gene regulatory network data as a CSV file.  
df.to_csv(x+"_rawRNA_celloracle_network_iGRN.csv")