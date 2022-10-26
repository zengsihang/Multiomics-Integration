import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import scglue
import seaborn as sns
from IPython import display
from matplotlib import rcParams
from networkx.algorithms.bipartite import biadjacency_matrix
from networkx.drawing.nx_agraph import graphviz_layout

scglue.plot.set_publication_params()
rcParams['figure.figsize'] = (4, 4)

rna = ad.read_h5ad("./data/rna-emb.h5ad")
atac = ad.read_h5ad("./data/atac-emb.h5ad")
guidance_hvf = nx.read_graphml("./data/guidance-hvf.graphml.gz")

rna.var["name"] = rna.var_names
atac.var["name"] = atac.var_names

genes = rna.var.query("highly_variable").index
peaks = atac.var.query("highly_variable").index

features = pd.Index(np.concatenate([rna.var_names, atac.var_names]))
feature_embeddings = np.concatenate([rna.varm["X_glue"], atac.varm["X_glue"]])

skeleton = guidance_hvf.edge_subgraph(
    e for e, attr in dict(guidance_hvf.edges).items()
    if attr["type"] == "fwd"
).copy()

reginf = scglue.genomics.regulatory_inference(
    features, feature_embeddings,
    skeleton=skeleton, random_state=0
)

gene2peak = reginf.edge_subgraph(
    e for e, attr in dict(reginf.edges).items()
    if attr["qval"] < 0.05
)

scglue.genomics.Bed(atac.var).write_bed("peaks.bed", ncols=3)
scglue.genomics.write_links(
    gene2peak,
    scglue.genomics.Bed(rna.var).strand_specific_start_site(),
    scglue.genomics.Bed(atac.var),
    "gene2peak.links", keep_attrs=["score"]
)

loc = rna.var.loc["Gad2"]
chrom = loc["chrom"]
chromLen = loc["chromEnd"] - loc["chromStart"]
chromStart = loc["chromStart"] - chromLen
chromEnd = loc["chromEnd"] + chromLen
print(chrom)
print(chromStart)
print(chromEnd)
# !pyGenomeTracks --tracks tracks.ini \
#     --region {chrom}:{chromStart}-{chromEnd} \
#     --outFileName tracks.png 2> .
# display.Image("tracks.png")