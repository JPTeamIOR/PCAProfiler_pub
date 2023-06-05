# Atlas generation scripts
Scripts to run in order to generate the Prostate Cancer Atlas.

### Requirements
The following requirements must be met:
 - Singularity ( Apptainer ) installed
 - Nextflow
 - fasterq-dump and prefetch ( sra-tools)
 - gdc-client
 - samtools
 - A folder (Utils directory) containing the tokens for restricted access databases under the paths ${BASE_UTILS}/tokens/dbGAP.ngc and ${BASE_UTILS}/tokens/TCGA.txt. The STAR reference will be generated in this folder with the first run and the references from GENCODE will be downloaded here.
 - In the folder where you run the 01_Download_and_process.sh there must be a subfolder called metadata with the information about the data to download, provided in this folder:
    - "${BASE_DIR}/metadata/gdc_manifest.TCGA_prostate.txt"
    - "${BASE_DIR}/metadata/SRR_ids.txt"
 - R with the requirements listed in the given scripts
## 01: Download and process

The first script download and process the raw FASTQ files downloaded from SRA or GDC ( TCGA ).
You can set the limits for each job at lines 52 and 53 with the variables 'threads' and 'max_mem'

### Run the script
To run all the samples, you can use the following commands:
```
./01_Download_and_process.sh SRA 0 100000 ~/Utils/
./01_Download_and_process.sh TCGA 0 100000 ~/Utils/
```
Where 0 is the lower limit, 100000 is the upper ( you can divide the work in batches ) and the ~/Utils/ is the BASE_UTILS folder containing the tokens.
The first process will generate a STAR reference in the Utils folder that will be reused by other processes, the first time wait that at least one of job is completed to run other batches.

### Output
A an ./output folder containing all the results. ( Salmon quantification, fastqc and trimgalore )

## 02: Collect data
Collect the data from the output folder and compress them. 

### Run the script
In the folder where you ran the previous script, and therefore containing the subfolder output, run
```
./02_Collect_data.sh
```

### Output
An export folder containing subfolders named with the sample ID and the compressed salmon's output quant.sf.gz.


# R scripts

Run the R scripts in order. Here's the sessionInfo() output:

```
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Rome
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] GSVA_1.48.0                 RColorBrewer_1.1-3          SingleCellExperiment_1.22.0 biomaRt_2.56.0              sva_3.48.0                  BiocParallel_1.34.0        
 [7] genefilter_1.82.0           mgcv_1.8-42                 nlme_3.1-162                DESeq2_1.40.0               SummarizedExperiment_1.30.0 Biobase_2.60.0             
[13] MatrixGenerics_1.12.0       matrixStats_0.63.0          GenomicRanges_1.52.0        GenomeInfoDb_1.36.0         IRanges_2.34.0              S4Vectors_0.38.0           
[19] BiocGenerics_0.46.0         tximport_1.28.0             edgeR_3.42.0                limma_3.56.0                jsonlite_1.8.4              ggrepel_0.9.3              
[25] lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0               dplyr_1.1.2                 purrr_1.0.1                 readr_2.1.4                
[31] tibble_3.2.1                tidyverse_2.0.0             tidyr_1.3.0                 ggplot2_3.4.2              

loaded via a namespace (and not attached):
 [1] DBI_1.1.3                 bitops_1.0-7              GSEABase_1.62.0           rlang_1.1.0               magrittr_2.0.3            compiler_4.3.0           
 [7] RSQLite_2.3.1             DelayedMatrixStats_1.22.0 png_0.1-8                 vctrs_0.6.2               pkgconfig_2.0.3           crayon_1.5.2             
[13] fastmap_1.1.1             dbplyr_2.3.2              XVector_0.40.0            labeling_0.4.2            utf8_1.2.3                tzdb_0.3.0               
[19] graph_1.78.0              bit_4.0.5                 beachmat_2.16.0           zlibbioc_1.46.0           cachem_1.0.7              slingshot_2.8.0          
[25] progress_1.2.2            blob_1.2.4                rhdf5filters_1.12.0       DelayedArray_0.25.0       Rhdf5lib_1.22.0           irlba_2.3.5.1            
[31] parallel_4.3.0            prettyunits_1.1.1         R6_2.5.1                  stringi_1.7.12            Rcpp_1.0.10               igraph_1.4.2             
[37] Matrix_1.5-1              splines_4.3.0             timechange_0.2.0          tidyselect_1.2.0          rstudioapi_0.14           codetools_0.2-19         
[43] curl_5.0.0                lattice_0.20-45           withr_2.5.0               KEGGREST_1.40.0           survival_3.5-5            BiocFileCache_2.8.0      
[49] TrajectoryUtils_1.8.0     xml2_1.3.4                Biostrings_2.68.0         BiocManager_1.30.20       pillar_1.9.0              filelock_1.0.2           
[55] generics_0.1.3            RCurl_1.98-1.12           hms_1.1.3                 sparseMatrixStats_1.12.0  munsell_0.5.0             scales_1.2.1             
[61] princurve_2.1.6           xtable_1.8-4              glue_1.6.2                tools_4.3.0               ScaledMatrix_1.7.1        annotate_1.78.0          
[67] locfit_1.5-9.7            XML_3.99-0.14             rhdf5_2.44.0              grid_4.3.0                AnnotationDbi_1.62.0      colorspace_2.1-0         
[73] GenomeInfoDbData_1.2.10   BiocSingular_1.16.0       HDF5Array_1.28.0          rsvd_1.0.5                cli_3.6.1                 rappdirs_0.3.3           
[79] fansi_1.0.4               gtable_0.3.3              digest_0.6.31             htmlwidgets_1.6.2         farver_2.1.1              memoise_2.0.1            
[85] htmltools_0.5.5           lifecycle_1.0.3           httr_1.4.5                bit64_4.0.5 

```
