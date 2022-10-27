# Multiomics-Integration
Course project of Computational Biology

## Download data
`wget http://download.gao-lab.org/GLUE/tutorial/Chen-2019-RNA.h5ad`

`wget http://download.gao-lab.org/GLUE/tutorial/Chen-2019-ATAC.h5ad`

`wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz`

## run codes


`python setup.py develop  # Build scglue`

`python data_preprocessing.py`

`python train.py`
