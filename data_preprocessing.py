import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
from matplotlib import rcParams
import matplotlib.pyplot as plt

rna = ad.read_h5ad("./data/Chen-2019-RNA.h5ad")
atac = ad.read_h5ad("./data/Chen-2019-ATAC.h5ad")
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
plt.close()

scglue.data.lsi(atac, n_components=100, n_iter=15)
sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
sc.tl.umap(atac)
sc.pl.umap(atac, color="cell_type")
plt.savefig("Chen-2019-ATAC-UMAP.png")
plt.close()

print(rna.var.head())
scglue.data.get_gene_annotation(
    rna, gtf="./data/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)
print(rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head())
print(rna.var.head())
print(rna.var.columns)

print(atac.var_names[:5])

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
print(atac.var.head())

guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
print(guidance)

scglue.graph.check_graph(guidance, [rna, atac])

print(atac.var.head())

rna.write("./data/rna-pp.h5ad", compression="gzip")
atac.write("./data/atac-pp.h5ad", compression="gzip")
nx.write_graphml(guidance, "./data/guidance.graphml.gz")

