import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams

rna = ad.read_h5ad("../dataset/Chen-2019-RNA.h5ad")
atac = ad.read_h5ad("../dataset/Chen-2019-ATAC.h5ad")
print(rna)
print(atac)

