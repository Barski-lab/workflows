##########################################################################################################
# Dockerfile
#
# Software:         set of R scripts for single-cell data analysis
# Software Version: v0.0.13
# Description:      Dockerized version of the set of R scripts for
#                   single-cell data analysis
# Website:          https://github.com/Barski-lab/workflows
# Provides:         set of R scripts for single-cell data analysis
# Base Image:       satijalab/seurat:4.0.6
# Build Cmd:        docker build --no-cache --rm -t biowardrobe2/sc-tools:v0.0.13 -f sc-tools-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/sc-tools:v0.0.13
# Run Cmd:          docker run --rm -ti biowardrobe2/sc-tools:v0.0.13 /bin/bash
##########################################################################################################
#
# Notes:
# - if you have problems compiling Harmony
#   https://github.com/immunogenomics/harmony/pull/165
#
# v0.0.13
# - Doesn't fail anymore when analyzing single dataset
# - Updated sc_rna_de_pseudobulk.R script to
#   * produce MDS plot
#   * to save HOPACH clustered aggregated GCT heatmap
#   * to use dittoHeatmap for cell heatmap
# - Updated sc_[rna/atac/wnn]_cluster.R scripts to allow
#   selecting cluster algorithm
# - Need to intall Harmony from GitHub as it disappeared
#   from BioConductor and official R repository
#
# v0.0.12
# - add sc_triangulate.R to run scTriangulate
# - install certain packages with manually provided
#   dependencies as R updates Seurat to the latest
#   version by default
# - need to have fixed versions for scipy and numpy
#   to be able to build numpy from source because of
#   seg. fault when running it from reticulate
#   https://stackoverflow.com/questions/70711946/reticulate-segfaults-with-call-to-plt-plot
#   https://github.com/rstudio/reticulate/issues/1133#issuecomment-1021783041
#
# v0.0.11
# - add integration with Harmony in sc_atac_reduce.R script
# - updated --barcodes parameter to extend metadata
# - updated --norm parameter in sc_atac_reduce.R
# - updated --callby parameter in sc_multiome_filter.R
#
# v0.0.10
# - add integration with Harmony on sc_rna_reduce.R script
# - bug fix with cell cycle removal in not integrated SCTransformed dataset
# - updated UCSC Cell Browser from 1.1.1 to 1.2.1
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
ENV SCTRIANGULATE_VERSION "0.12.0"
ENV SCIPY_VERSION "1.7.3"
ENV NUMPY_VERSION "1.22.0"

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
COPY ./scripts/sc_tools/sc_triangulate.R /usr/local/bin/sc_triangulate.R
COPY ./scripts/sc_tools/modules/*.R /usr/local/bin/modules/

### Installing dependencies
RUN apt-get update && \
    apt-get install libgcc-10-dev python3-dev python3-pip libxml2-dev libcurl4-openssl-dev \
                    libssl-dev pandoc libudunits2-dev libgdal-dev libcairo2-dev libharfbuzz-dev \
                    libfribidi-dev libbz2-dev tk libfontconfig1-dev libfreetype6-dev libpng-dev \
                    libtiff5-dev libjpeg-dev vim -y && \
    pip3 install --no-binary="numpy" numpy==${NUMPY_VERSION} --ignore-installed && \
    pip3 install scipy==${SCIPY_VERSION} && \
    pip3 install anndata && \
    pip3 install cellbrowser==${CB_VERSION} && \
    pip3 install -U numba && \
    pip3 install MACS2==${MACS2_VERSION} && \
    pip3 install sctriangulate==${SCTRIANGULATE_VERSION} && \
    R -e 'install.packages("devtools", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("argparse", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("tidyverse", repo = "https://cloud.r-project.org/")' && \
    R -e 'install.packages("modules", repo = "https://cloud.r-project.org/")' && \
    R -e "BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', \
                                 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', \
                                 'batchelor', 'Matrix.utils', 'DESeq2', 'Rsamtools', 'rtracklayer', \
                                 'glmGamPoi', 'Nebulosa', 'limma', 'EnhancedVolcano', 'LoomExperiment', \
                                 'hopach', 'Glimma', 'dittoSeq', 'cmapR'))" && \
    R -e 'install.packages("fastmatch")' && \
    R -e 'install.packages("RcppRoll")' && \
    R -e 'install.packages("Signac", repo = "https://cloud.r-project.org/", dependencies = FALSE)' && \
    R -e 'devtools::install_github("cellgeni/sceasy", upgrade = "never")' && \
    R -e 'devtools::install_github("KlugerLab/DAseq", upgrade = "never")' && \
    R -e 'install.packages("LambertW")' && \
    R -e 'install.packages("nortest")' && \
    R -e 'install.packages("doParallel")' && \
    R -e 'install.packages("doRNG")' && \
    R -e 'install.packages("butcher")' && \
    R -e 'install.packages("recipes")' && \
    R -e 'install.packages("bestNormalize", repo = "https://cloud.r-project.org/", dependencies = FALSE)' && \
    R -e 'devtools::install_github("immunogenomics/harmony", upgrade = "never")' && \
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
    chmod +x /usr/local/bin/sc_triangulate.R && \
### Cleaning
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true
