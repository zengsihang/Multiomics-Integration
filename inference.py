from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
import matplotlib.pyplot as plt

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

rna = ad.read_h5ad("data/rna-pp.h5ad")
atac = ad.read_h5ad("data/atac-pp.h5ad")
guidance = nx.read_graphml("data/guidance.graphml.gz")
scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)
scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)
guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

glue = scglue.models.load_model("glue.dill")
dx = scglue.models.integration_consistency(
    glue, {"rna": rna, "atac": atac}, guidance_hvf
)
print(dx)

_ = sns.lineplot(x="n_meta", y="consistency", data=dx).axhline(y=0.05, c="darkred", ls="--")
plt.savefig("Chen-2019-SCGLUE-Integration-Consistency.png")
plt.close()

rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)

combined = ad.concat([rna, atac])
sc.pp.neighbors(combined, use_rep="X_glue", metric="cosine")
sc.tl.umap(combined)
sc.pl.umap(combined, color=["cell_type", "domain"], wspace=0.65)

feature_embeddings = glue.encode_graph(guidance_hvf)
feature_embeddings = pd.DataFrame(feature_embeddings, index=glue.vertices)
print(feature_embeddings.iloc[:5, :5])

rna.varm["X_glue"] = feature_embeddings.reindex(rna.var_names).to_numpy()
atac.varm["X_glue"] = feature_embeddings.reindex(atac.var_names).to_numpy()

rna.write("data/rna-emb.h5ad", compression="gzip")
atac.write("data/atac-emb.h5ad", compression="gzip")
nx.write_graphml(guidance_hvf, "data/guidance-hvf.graphml.gz")