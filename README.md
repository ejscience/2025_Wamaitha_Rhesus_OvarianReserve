# 2024_Wamaitha_Rhesus_OvarianReserve
Github for Wamaitha, S. article: "Defining the cell and molecular origins of the ovarian reserve" 

## Overview

This project is focused on running a full mixed probe analysis, as detailed in the `20231202_MixedProbesFullRun.nb.html` file. 
`20231202_MixedProbesFullRun` : The HTML document contains the full experimental run and includes the following sections

## File Structure
2024_Wamaitha_Rhesus_OvarianReserve/
├── README.md
├── spaceranger_array_submission.sh
├── 20231202_MixedProbesFullRun.nb.html
├── UpdatedSpaceRanger/
│   ├── MakingMixedProbes.sh
│   ├── MakingSingleProbes.R
│   ├── singlehit_mixed_single.probes.csv
│   ├── single.probe.csv
│   └── RhesusMacaqueProbeLists/
│       ├── FRP_Human_probes_on_Macaque_1.csv
│       ├── FRP_Human_probes_on_Macaque_2.csv
│       ├── FRP_Human_probes_on_Macaque_3.csv
│       └── FRP_Human_probes_on_Macaque_4.csv

- Folder `UpdatedSpaceRanger` contains the files below which will make the probe set list to align the data.
    - `UpdatedSpaceRanger/MakingSingleProbes.R`: This is the file that will be run first to make the single hits of the Rhesus Mulatta using human probes
    - `UpdatedSpaceRanger/MakingMixedProbes.sh`: This is the file that will then take the single hits and add in the mixed probes
    - The rest of the `.csv` files are the probes that were added from running the above files. And the folder `UpdatedSpaceRanger/RhesusMacaqueProbeLists` contains the 10X genomics hits that they found.
- `spaceranger_array_submission.sh`: This is the script to run spaceranger on a cluster but clearly shows how to integrate the new probe set.
- `20231202_MixedProbesFullRun.nb.html`: The main output of the analysis, generated using the notebook `20231202_MixedProbesFullRun.Rmd`. This file contains the detailed steps, code, and results from the experiment.

## Usage

1. **Viewing Results**: Open the `20231202_MixedProbesFullRun.nb.html` file in any web browser to view the interactive analysis, which includes plots and tables.
2. **Reproducing Analysis**:
    - Clone the repository.
    - Start by running the file `UpdatedSpaceRanger/MakingSingleProbes.R` and `UpdatedSpaceRanger/MakingMixedProbes.sh`. These will make the probe lists that will be used in spaceranger
    - Then follow the steps in our `spaceranger_array_submission.sh` to include the correct probe list to your spaceranger submission.
    - Take your 
    - Follow the code snippets in the `Methods` section to set up the experimental environment and run your own version of the analysis.
    - Ensure all required libraries and dependencies are installed, these are at the end in the SessionInfo() of the HTML `20231202_MixedProbesFullRun.nb.html` or see below.

## Requirements

To rerun the analysis, make sure the following dependencies are installed:
```
sessionInfo()
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.2.0             harmony_1.2.0               Rcpp_1.0.13                 DoubletFinder_2.0.4        
 [5] SoupX_1.6.2                 BiocManager_1.30.24         readxl_1.4.3                clustree_0.5.1             
 [9] ggraph_2.2.1                RCurl_1.98-1.16             cowplot_1.1.3               Matrix_1.6-5               
[13] ggsignif_0.6.4              data.table_1.15.4           reshape_0.8.9               biomaRt_2.58.2             
[17] scales_1.3.0                org.Hs.eg.db_3.18.0         AnnotationDbi_1.64.1        clusterProfiler_4.10.1     
[21] RColorBrewer_1.1-3          pheatmap_1.0.12             ggdendro_0.2.0              ggrepel_0.9.5              
[25] DESeq2_1.42.1               SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0      
[29] matrixStats_1.3.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0             
[33] S4Vectors_0.40.2            BiocGenerics_0.48.1         ggpattern_1.1.1             ComplexHeatmap_2.18.0      
[37] lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
[41] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
[45] ggplot2_3.5.1               tidyverse_2.0.0             colorRamp2_0.1.0            ggprism_1.0.5              
[49] scCustomize_2.1.2           Seurat_5.1.0                SeuratObject_5.0.2          sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] progress_1.2.3            goftest_1.2-3             Biostrings_2.70.3         vctrs_0.6.5               spatstat.random_3.3-1    
  [6] digest_0.6.37             png_0.1-8                 shape_1.4.6.1             deldir_2.0-4              parallelly_1.38.0        
 [11] MASS_7.3-59               reshape2_1.4.4            httpuv_1.6.15             foreach_1.5.2             qvalue_2.34.0            
 [16] withr_3.0.1               ggrastr_1.0.2             xfun_0.47                 ggfun_0.1.5               survival_3.7-0           
 [21] memoise_2.0.1             ggbeeswarm_0.7.2          janitor_2.2.0             gson_0.1.0                tidytree_0.4.6           
 [26] zoo_1.8-12                GlobalOptions_0.1.2       pbapply_1.7-2             prettyunits_1.2.0         rematch2_2.1.2           
 [31] KEGGREST_1.42.0           promises_1.3.0            httr_1.4.7                globals_0.16.3            fitdistrplus_1.2-1       
 [36] rstudioapi_0.16.0         miniUI_0.1.1.1            generics_0.1.3            DOSE_3.28.2               curl_5.2.2               
 [41] zlibbioc_1.48.2           polyclip_1.10-7           GenomeInfoDbData_1.2.11   SparseArray_1.2.4         xtable_1.8-4             
 [46] doParallel_1.0.17         evaluate_0.24.0           S4Arrays_1.2.1            BiocFileCache_2.10.2      hms_1.1.3                
 [51] irlba_2.3.5.1             colorspace_2.1-1          filelock_1.0.3            hdf5r_1.3.11              ROCR_1.0-11              
 [56] reticulate_1.38.0         spatstat.data_3.1-2       magrittr_2.0.3            lmtest_0.9-40             snakecase_0.11.1         
 [61] later_1.3.2               viridis_0.6.5             ggtree_3.10.1             lattice_0.22-6            glmGamPoi_1.14.3         
 [66] spatstat.geom_3.3-2       future.apply_1.11.2       scattermore_1.2           XML_3.99-0.17             shadowtext_0.1.4         
 [71] RcppAnnoy_0.0.22          pillar_1.9.0              nlme_3.1-166              iterators_1.0.14          compiler_4.3.0           
 [76] RSpectra_0.16-2           stringi_1.8.4             tensor_1.5                plyr_1.8.9                crayon_1.5.3             
 [81] abind_1.4-5               gridGraphics_0.5-1        locfit_1.5-9.10           graphlayouts_1.1.1        bit_4.0.5                
 [86] fastmatch_1.1-4           codetools_0.2-20          bslib_0.8.0               paletteer_1.6.0           GetoptLong_1.0.5         
 [91] plotly_4.10.4             mime_0.12                 splines_4.3.0             circlize_0.4.16           fastDummies_1.7.4        
 [96] dbplyr_2.5.0              sparseMatrixStats_1.14.0  prismatic_1.1.2           HDO.db_0.99.1             cellranger_1.1.0         
[101] knitr_1.48                blob_1.2.4                utf8_1.2.4                clue_0.3-65               fs_1.6.4                 
[106] listenv_0.9.1             DelayedMatrixStats_1.24.0 ggplotify_0.1.2           tzdb_0.4.0                tweenr_2.0.3             
[111] pkgconfig_2.0.3           tools_4.3.0               cachem_1.1.0              RSQLite_2.3.7             viridisLite_0.4.2        
[116] DBI_1.2.3                 fastmap_1.2.0             rmarkdown_2.28            ica_1.0-3                 sass_0.4.9               
[121] dotCall64_1.1-1           RANN_2.6.2                farver_2.1.2              tidygraph_1.3.1           scatterpie_0.2.3         
[126] yaml_2.3.10               cli_3.6.3                 leiden_0.4.3.1            lifecycle_1.0.4           uwot_0.2.2               
[131] BiocParallel_1.36.0       timechange_0.3.0          gtable_0.3.5              rjson_0.2.21              ggridges_0.5.6           
[136] progressr_0.14.0          parallel_4.3.0            ape_5.8                   jsonlite_1.8.8            RcppHNSW_0.6.0           
[141] bitops_1.0-8              bit64_4.0.5               Rtsne_0.17                yulab.utils_0.1.7         spatstat.utils_3.1-0     
[146] jquerylib_0.1.4           GOSemSim_2.28.1           spatstat.univar_3.0-0     lazyeval_0.2.2            shiny_1.9.1              
[151] htmltools_0.5.8.1         enrichplot_1.22.0         GO.db_3.18.0              sctransform_0.4.1         rappdirs_0.3.3           
[156] glue_1.7.0                spam_2.10-0               XVector_0.42.0            treeio_1.26.0             gridExtra_2.3            
[161] igraph_2.0.3              R6_2.5.1                  labeling_0.4.3            cluster_2.1.6             aplot_0.2.3              
[166] DelayedArray_0.28.0       tidyselect_1.2.1          vipor_0.4.7               ggforce_0.4.2             xml2_1.3.6               
[171] future_1.34.0             munsell_0.5.1             KernSmooth_2.23-24        htmlwidgets_1.6.4         fgsea_1.28.0             
[176] rlang_1.1.4               spatstat.sparse_3.1-0     spatstat.explore_3.3-2    fansi_1.0.6               Cairo_1.6-2              
[181] beeswarm_0.4.0   
```

## Contact

For any questions or further information, feel free to reach out to the project lead.

---

*Updated on: October 21 2024*
