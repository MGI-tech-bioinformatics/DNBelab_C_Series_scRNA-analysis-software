#!env python3
import sys
import numpy as np
import pandas as pd
import scanpy as sc

### input data
data = sc.read(sys.argv[1],cache=False)
result_dir = sys.argv[2]
adata = data.T
adata.var_names_make_unique()
### filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
### for each cell compute fraction of counts in mito genes vs all genes
### add the total counts per cell as observations-annotation to adata
mito_genes = adata.var_names.str.startswith('mt-')
adata.obs['percent_mi'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
adata.obs['n_counts'] = adata.X.sum(axis=1)
############## Plot ###################
#sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)
#sc.pl.scatter(adata, x='n_counts', y='percent_mito')
#sc.pl.scatter(adata, x='n_counts', y='n_genes')

### real filter
adata = adata[adata.obs['n_genes'] > 200, :]
#adata = adata[adata.obs['n_genes'] < 4000, :]
adata = adata[adata.obs['percent_mi'] < 0.5, :]
### normalize
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata
with open("%s/cell_report.csv" % result_dir, "a") as o:
    o.write("Number of cells used for clustering,%s\nTotal Genes Detected,%s" % (adata.n_obs, adata.n_vars))
o.close()

###差异表达基因
sc.pp.highly_variable_genes(adata)
#sc.pl.highly_variable_genes(adata)
adata = adata[:, adata.var['highly_variable']]

sc.pp.regress_out(adata, ['n_counts', 'percent_mi'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.louvain(adata)
sc.tl.umap(adata)
#sc.tl.tsne(adata)

cluster=pd.concat([adata.obs['louvain'][:],adata.obs['percent_mi'][:],adata.obs['n_counts'][:],adata.obs['n_genes'][:]],axis=1)
cluster["UMAP_1"]=adata.obsm["X_umap"][:,0]
cluster["UMAP_2"]=adata.obsm["X_umap"][:,1]
cluster.columns = ["Cluster","percent.mi","nUMI","nGene","UMAP_1","UMAP_2"]
cluster = cluster[["Cluster","UMAP_1","UMAP_2","percent.mi","nUMI","nGene"]]
cluster.to_csv('%s/cluster.csv' % result_dir)

sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon', n_genes=50)
sc.tl.filter_rank_genes_groups(adata, min_fold_change=0.05, min_in_group_fraction=0.05)
df=pd.DataFrame(np.arange(len(adata.uns['rank_genes_groups_filtered']['scores'][0])*50*6).reshape((len(adata.uns['rank_genes_groups_filtered']['scores'][0])*50,6)),columns=['genes','cluster','p_val','avg_logFC','scores','p_val_adj'],index=range(1,len(adata.uns['rank_genes_groups_filtered']['scores'][0])*50+1))
for clu in range(len(adata.uns['rank_genes_groups_filtered']['scores'][0])):
    for row in range(50):
        df.loc[clu*50+row+1]=[adata.uns['rank_genes_groups_filtered']['names'][row][clu],clu,adata.uns['rank_genes_groups_filtered']['pvals'][row][clu],adata.uns['rank_genes_groups_filtered']['logfoldchanges'][row][clu],adata.uns['rank_genes_groups_filtered']['scores'][row][clu],adata.uns['rank_genes_groups_filtered']['pvals_adj'][row][clu]]
df=df[~df['genes'].isin([np.nan])]
df.index = range(1,len(df)+1)
df.to_csv('%s/marker.csv' % result_dir)
adata.write("%s/cluster.h5ad" % result_dir)

##                                    louvain  percent_mito  n_counts  n_genes
##TCAGTACACGAGCAAGCGAC-V300036264_L01       9      0.006178    7932.0     3784
##ATATGATGCTTATTACTCTC-V300036264_L01       9      0.010308    6791.0     3307
