import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gprofiler import GProfiler

# Load dataset (modify with actual dataset path or GEO dataset)
adata = sc.read_10x_mtx("path_to_data")  # Replace with dataset path

# Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]

# PCA and clustering
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# Identify marker genes for T and B cells
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
marker_genes = ['CD3D', 'CD3E', 'CD19', 'MS4A1']  # T and B cell markers
sc.pl.umap(adata, color=marker_genes)

# Functional Enrichment Analysis (GO, KEGG)
gp = GProfiler(return_dataframe=True)
enrich_res = gp.profile(organism='hsapiens', query=marker_genes)

# Plot Enrichment Results
sns.barplot(data=enrich_res, x='p_value', y='name')
plt.xlabel("P-value")
plt.ylabel("Pathway")
plt.title("Functional Enrichment Analysis")
plt.show()

# Save results
enrich_res.to_csv("enrichment_results.csv")
adata.write("processed_scRNA.h5ad")
