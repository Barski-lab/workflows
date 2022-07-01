##########################################################################################################
# Dockerfile
#
# Software:         set of R scripts for single-cell data analysis
# Software Version: v0.0.9
# Description:      Dockerized version of the set of R scripts for
#                   single-cell data analysis
# Website:          https://github.com/Barski-lab/workflows
# Provides:         set of R scripts for single-cell data analysis
# Base Image:       satijalab/seurat:4.0.6
# Build Cmd:        docker build --no-cache --rm -t biowardrobe2/sc-tools:v0.0.9 -f sc-tools-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/sc-tools:v0.0.9
# Run Cmd:          docker run --rm -ti biowardrobe2/sc-tools:v0.0.9 /bin/bash
##########################################################################################################
#
# v0.0.9
# - export h5ad files
# - minor graphics correction
# - add sc_rna_da_cells.R script
#
# v0.0.8
# - rename sc_de_pseudobulk.R to sc_rna_de_pseudobulk.R
#   added more UMAPs to plot
#
# v0.0.7
# - added sc_de_pseudobulk.R script for pseudobulk DE analysis
# - install limma package
#
# v0.0.6
# - added sc_ctype_assign.R script for manual cell types
#   assignment
#
# v0.0.5
# - updated sc_[rna/wnn]_cluster.R script to use densoty plots
#   from Nebulosa
# - removed --regressrnaumi from sc_rna_reduce.R
# - made --regressgenesin sc_rna_reduce.R to use list of genes
#
# v0.0.4
# - add/update scripts
#   sc_rna_filter.R
#   sc_multiome_filter.R
#   sc_rna_reduce.R
#   sc_atac_reduce.R
#   sc_rna_cluster.R
#   sc_atac_cluster.R
#   sc_wnn_cluster.R
# - install MACS2 after numba as it fails on numpy
#
# v0.0.3
# - add sc_gex_reduce.R script
# - add sc_gex_cluster.R
# Notes:
# - pip3 install -U numba is required. This will reinstall numpy
#   to the version that is needed by numba thus making possible to
#   use umap in R through reticulate (use umap-learn in runUMAP).
#
# v0.0.2
# - add sc_gex_filter.R script
#
# v0.0.1
# - initial version of sc_multiome_filter.R script
#
#################################################################


### Base Image
FROM satijalab/seurat:4.0.6
LABEL maintainer="misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive


################## BEGIN INSTALLATION ######################

WORKDIR /tmp

ENV R_MAX_VSIZE=200000000000
ENV CB_VERSION "1.1.1"
ENV MACS2_VERSION "2.2.7.1"

COPY ./scripts/sc_tools/sc_rna_filter.R /usr/local/bin/sc_rna_filter.R
COPY ./scripts/sc_tools/sc_multiome_filter.R /usr/local/bin/sc_multiome_filter.R
COPY ./scripts/sc_tools/sc_rna_reduce.R /usr/local/bin/sc_rna_reduce.R
COPY ./scripts/sc_tools/sc_atac_reduce.R /usr/local/bin/sc_atac_reduce.R
COPY ./scripts/sc_tools/sc_rna_cluster.R /usr/local/bin/sc_rna_cluster.R
COPY ./scripts/sc_tools/sc_atac_cluster.R /usr/local/bin/sc_atac_cluster.R
COPY ./scripts/sc_tools/sc_wnn_cluster.R /usr/local/bin/sc_wnn_cluster.R
COPY ./scripts/sc_tools/sc_ctype_assign.R /usr/local/bin/sc_ctype_assign.R
COPY ./scripts/sc_tools/sc_rna_de_pseudobulk.R /usr/local/bin/sc_rna_de_pseudobulk.R
COPY ./scripts/sc_tools/sc_rna_da_cells.R /usr/local/bin/sc_rna_da_cells.R
COPY ./scripts/sc_tools/modules/*.R /usr/local/bin/modules/

### Installing dependencies
RUN apt-get update && \
    apt-get install libgcc-10-dev python3-dev python3-pip libxml2-dev libcurl4-openssl-dev \
                    libssl-dev pandoc libudunits2-dev libgdal-dev libcairo2-dev libharfbuzz-dev \
                    libfribidi-dev libbz2-dev tk -y && \
    pip3 install scipy numpy anndata && \
    pip3 install cellbrowser==${CB_VERSION} && \
    pip3 install -U numba && \
    pip3 install MACS2==${MACS2_VERSION} && \
    R -e 'install.packages("devtools", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("argparse", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("tidyverse", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("modules", repo = "https://cloud.r-project.org/")' && \
    R -e "BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'DESeq2', 'Rsamtools', 'rtracklayer', 'glmGamPoi', 'Nebulosa', 'limma', 'EnhancedVolcano', 'LoomExperiment'))" && \
    R -e 'install.packages("Signac", repo = "https://cloud.r-project.org/")' && \
    R -e 'devtools::install_github("cellgeni/sceasy")' && \
    R -e 'devtools::install_github("KlugerLab/DAseq")' && \
    R -e 'install.packages("bestNormalize", repo = "https://cloud.r-project.org/")' && \
### Installing R scripts
    chmod +x /usr/local/bin/sc_rna_filter.R && \
    chmod +x /usr/local/bin/sc_multiome_filter.R && \
    chmod +x /usr/local/bin/sc_rna_reduce.R && \
    chmod +x /usr/local/bin/sc_atac_reduce.R && \
    chmod +x /usr/local/bin/sc_rna_cluster.R && \
    chmod +x /usr/local/bin/sc_atac_cluster.R && \
    chmod +x /usr/local/bin/sc_wnn_cluster.R && \
    chmod +x /usr/local/bin/sc_ctype_assign.R && \
    chmod +x /usr/local/bin/sc_rna_de_pseudobulk.R && \
    chmod +x /usr/local/bin/sc_rna_da_cells.R && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
