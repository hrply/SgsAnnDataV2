# How to install SgsAnnData in conda

The correct steps to install SgsAnnData without repeating bugs of dependencies: 

## Install R-base=4.4.1~4.4.3 and essential packages

It's recommended to install R-base 4.4.1 and then upgrade to 4.4.3 as the renv.lock of ArchR is builded in R 4.4.1.

```
conda install -c conda-forge r-base=4.4.3 radian rhash r-essentials compilers make gcc cmake wheel curl
conda install -c conda-forge xz bzip2 zlib xorg-x11-proto-devel-cos7-x86_64 libx11-devel-cos7-x86_64 libx11-cos7-x86_64 libxext-devel-cos7-x86_64 libxext-cos7-x86_64 cairo-devel-cos7-x86_64 cairo-cos7-x86_64 libxt-devel-cos7-x86_64 libzlib  libiconv libxml2-devel-cos7-x86_64 glpk r-glpkapi r-rglpk libgit2 r-git2r gsl r-gsl udunits2 libudunits2 r-udunits2 pandoc gdal libgdal libgdal-core-devel imagemagick r-devtools r-remotes
```

## Install ArchR and Seurat with environment packages by renv
### install renv and environment packages

Create renv project file and run R in this path

```
conda create -n sgs  # Assuming the conda environment name of SgsAnnData is sgs
conda activate sgs
mkdir -p $CONDA_PREFIX/renv && cd $CONDA_PREFIX/renv 
R
```

Then install renv and packages

```
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes") #if not install r-remotes
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools") #if not install r-devtools
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.20")
install.packages("renv")
library(renv)
setwd("$CONDA_PREFIX/renv")
download.file(url = "https://pub-9ae435458ecc412abbbc9420a502ec38.r2.dev/renv.lock", destfile = "./renv.lock")
renv::init(settings = list(sandbox.enabled = FALSE))
renv::settings$external.libraries("$CONDA_PREFIX/lib/R/library")
renv::init()
```

Tips: If some packages install failed, generally, errors occurs when "cannot find xx.so or xx.h", so you can read the error messages and install the nessary library packages. 

### install ArchR and Seurat

Then install ArchR & Seurat in R:

```
install.packages('Seurat')
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
q()
```

## Install Giotto

For Linux, there are several prerequisite installs: GDAL (>= 2.2.3), GEOS (>= 3.4.0), PROJ (>= 4.9.3), sqlite3

```
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo apt-get update
sudo apt-get install libgdal-dev libgeos-dev libproj-dev libtbb-dev libnetcdf-dev
conda install -c conda-forge gdal libgdal geos proj sqlite libsqlite
```

Then run in R

```
install.packages("terra") # if terra version is lower than 1.8-21, please run ``install.packages("https://cloud.r-project.org/src/contrib/Archive/terra/terra_1.8-29.tar.gz", repos = NULL, type = "source")
if(!"pak" %in% installed.packages()) install.packages("pak")
pak::pkg_install("drieslab/Giotto")
library("Giotto") # and install packages according to returned messages
```

## Install Signac

Run in R:

```
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac") # or "conda install -c bioconda r-signac"
```

**Installing genome assembly and gene annotation packages if nessary**

```
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
```
## Install Seurat supporting packages if nessary

Run in R

```
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
remotes::install_github("satijalab/seurat-data", force = TRUE)
remotes::install_github("satijalab/seurat-wrappers", force = TRUE)
remotes::install_github("satijalab/azimuth", ref = 'master')
```

## SgsAnnData installation

Run in R

```
BiocMamager::install(c("imager", "magick") # This two packages are easy to return error, it's recommended to install them separately.
devtools::install_github("bio-xtt/SgsAnnDataV2")
```

> ### for github packages install

> Use the code to add token to R
> ```
> install.packages("gitcreds")
> usethis::create_github_token()
> gitcreds::gitcreds_set()
> ```

**Remember install all packages in renv environment!** 

After finishing installation of all packages, update the lock file.

```
all_packages <- installed.packages()[, "Package"]
renv::snapshot(packages = all_packages)
```

## Experimental

Updating the BiocManager to 3.20:

1. Open the lockfile and change the repos like:
>Before
```
{
  "R": {
    "Version": "4.4.3",
    "Repositories": [
      {
        "Name": "BioCsoft",
        "URL": "https://bioconductor.org/packages/3.19/bioc"
      },
      {
        "Name": "BioCann",
        "URL": "https://bioconductor.org/packages/3.19/data/annotation"
      },
      {
        "Name": "BioCexp",
        "URL": "https://bioconductor.org/packages/3.19/data/experiment"
      },
      {
        "Name": "BioCworkflows",
        "URL": "https://bioconductor.org/packages/3.19/workflows"
      },
      {
        "Name": "BioCbooks",
        "URL": "https://bioconductor.org/packages/3.19/books"
      },
      {
        "Name": "CRAN",
        "URL": "https://p3m.dev/cran/2024-10-30"
      }
    ]
  },
  "Bioconductor": {
    "Version": "3.19"
  },
```
>After
```
{
  "R": {
    "Version": "4.4.3",
    "Repositories": [
      {
        "Name": "BioCsoft",
        "URL": "https://bioconductor.org/packages/3.20/bioc"
      },
      {
        "Name": "BioCann",
        "URL": "https://bioconductor.org/packages/3.20/data/annotation"
      },
      {
        "Name": "BioCexp",
        "URL": "https://bioconductor.org/packages/3.20/data/experiment"
      },
      {
        "Name": "BioCworkflows",
        "URL": "https://bioconductor.org/packages/3.20/workflows"
      },
      {
        "Name": "BioCbooks",
        "URL": "https://bioconductor.org/packages/3.20/books"
      },
      {
        "Name": "CRAN",
        "URL": "https://mirrors.zju.edu.cn/CRAN" # change to your prefered repos
      }
    ]
  },
  "Bioconductor": {
    "Version": "3.20"
  },
```

2. restart R and delete the old version.
    ``renv::remove("BiocManager")``
3. install the new version. 
    ```
    renv::install("BiocManager@1.3.25")
    BiocManager::install(version = "3.20")
    ```
4. update other packages if nessary. ``BiocManager::install()``

---
---

# forked from bio-xtt
# SgsAnnData
SgsAnnData is an R package that facilitates the seamless conversion of single-cell analysis object from popular tools such as Seurat, Giotto, Signac, and ArchR into the AnnData. This format can be directly visualized in SGS which is an interactive browser for single-cell and spatial multimodal datasets.

![](https://211.bioinfotoolkits.net:10290/sgs/website/feature_4.png)

> ## Warning
> Please note that the usage of SGS requires the prior installation of R packages such as Seurat, Giotto, Signac, and anndata. It is important to ensure that these packages are installed beforehand for smooth execution of SGS workflows.

## Installation
Currently, SgsAnndata can only be installed through GitHub. However, we are actively working on providing a Bioconductor installation method in the near future. 
```
devtools::install_github("bio-xtt/SgsAnnDataV2")
```
## Run SgsAnnData
```
renv::activate()
library(anndata)
```

## Usage
### Seurat to AnnData
### Single-Cell Transcriptome Object Conversion
When Seurat contains multiple assays, users can provide multi assay names they wish to export. This will automatically generates the h5ad for each assay type, like RNA.h5ad,integration.h5ad etc. We also offer the flexibility for users to output marker dataframes as a list, with the corresponding information stored in **adata.uns['markers']**.
```
library(Seurat)
library(anndata)
SeuratToAnndata(object=scRNA,
                outpath="/test_adata",
                assays=c("RNA", "SCT"),
                groups = NULL,
                reductions = c("tsne","umap"),
                markersDF = list("RNA"=rna_marker_df, "SCT"=sct_marker_df)) 
```

### Spatial Transcriptome Object Conversion
When converting spatial transcriptomic data, we automatically determine the spatial data type based on the object's structure. We return the corresponding **h5ad** file based on the identified type. The spatial coordinates are stored within **adata.obsm['spatial']**. If the analyzed object includes spatially organized slice information, we store it in **adata.uns['spatial']**. Multiple slices are differentiated using different **"library_id"**, ensuring clear distinction between them.
```
SeuratToAnndata(object=scRNA,
                outpath="/test_adata",
                assays=c("RNA", "Spatial"),
                groups = c("seurat_clusters", "region"),
                reductions = "tsne",
                markersDF = list("RNA"=rna_marker_df, "Spatial"=st_marker_df)) 
```

### Signac to Anndata
We recommend exporting **peak matrix**, **genescore matrix**, and **motif matrix** from single-cell ATAC data analysis as h5ad files. If desired, these files can be further converted into [**Mudata**](https://mudata.readthedocs.io/en/latest/), which offers a comprehensive and annotated multimodal dataset structure. Additionally, users can use **export_pwm** to output results from co-accessibility analysis and motif enrichment analysis, respectively. Users can access motif-related information via **adata.uns.["motifs"]**.
```
SignacToAnndata(object=scATAC,
                outpath = "/test_adata",
                assays=c("RNA", "Peaks", "Motif"),
                assay.types=("gene","peak","motif"),
                markersDF = list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
                groups = c("clusters","sample", "age"),
                reductions = c("umap","tsne"),
                export_pwm = TRUE)   # to export motif pwm 
```

### ArchR to AnnData
Similar to the Signac object conversion, users can convert ArchR objects and export relevant information based on their needs.
```
ArchrToAnndata(object=project5,
               outpath="/test_adata",
               assays=c("RNA", "Peaks", "Motif"),
               assay.types=("gene","peak","motif"),
               markersDF=list("RNA"=gene_marker_df, "Peaks"=marker_peak_df, "Motif"=marker_motif_df),
               groups = NULL, ## export all cell meta
               reductions = NULL, ## export all reduction coords
               export_pwm =TRUE) 
```

### Giotto to AnnData
Similar to Seurat's conversion of spatial transcriptomic analysis results, we seamlessly incorporate tissue slices and spatial coordinates from Giotto object into the output h5ad object by default.
```
giottoToAnnData(object = giotto,
                outpath = "/test_adata",
                markerDF = NULL)
```

## Citiation
Xia, T., Sun, J., Lu, F., Luo, Y., Mao, Y., Xu, L., & Wang, Y. (2024). Empowering Integrative and Collaborative Exploration of Single-Cell and Spatial Multimodal Data with SGS. bioRxiv, (), 2024.07.19.604227. Accessed July 23, 2024. https://doi.org/10.1101/2024.07.19.604227.










