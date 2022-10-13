import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import matplotlib.pyplot as plt

rna = ad.read_h5ad("../dataset/Chen-2019-RNA.h5ad")
atac = ad.read_h5ad("../dataset/Chen-2019-ATAC.h5ad")
print(rna)
print(atac)

rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")
sc.pp.neighbors(rna, metric="cosine")
sc.tl.umap(rna)
sc.pl.umap(rna, color="cell_type")
plt.savefig("Chen-2019-RNA-UMAP.png")