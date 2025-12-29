# Hemizygosity and HR deficiency drive resistance in breast cancer

## Safonov *et al.*
### https://someurl
&nbsp;
&nbsp;
&nbsp;

![Front page](https://github.com/antonmsafonov/germline-somatic-RB1/blob/repack/etc/splash.png)

#### Clone repository
```
git clone https://github.com/antonmsafonov/germline-somatic-RB1.git
cd germline-somatic-RB1/R
```

#### Render markdowns in R
```
library("rmarkdown")
rmarkdown::render(input = "Figure_1.Rmd", output_dir = "../res/")
rmarkdown::render(input = "Figure_2.Rmd", output_dir = "../res/")
rmarkdown::render(input = "Figure_3.Rmd", output_dir = "../res/")
rmarkdown::render(input = "Figure_4.Rmd", output_dir = "../res/")
```

#### R session info
```
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS  10.16

Matrix products: default
BLAS/LAPACK: /Users/davidbrown/.miniconda3/lib/libopenblasp-r0.3.10.dylib

locale:
[1] C

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] viridis_0.5.1                     viridisLite_0.3.0                
 [3] forcats_0.5.0                     stringr_1.4.0                    
 [5] purrr_0.3.4                       tibble_3.0.3                     
 [7] tidyverse_1.3.0                   tidyr_1.1.1                      
 [9] survminer_0.4.9                   ggpubr_0.2                       
[11] survival_3.2-10                   superheat_1.0.0                  
[13] reshape2_1.4.4                    readr_1.3.1                      
[15] RColorBrewer_1.1-2                pander_0.6.3                     
[17] mime_0.9                          magrittr_1.5                     
[19] maftools_2.12.05                  lubridate_1.7.9                  
[21] logistf_1.26.1                    knitr_1.29                       
[23] jsonlite_1.7.0                    htmltools_0.5.0                  
[25] Gviz_1.32.0                       gridExtra_2.3                    
[27] ggsignif_0.6.0                    ggrepel_0.8.2                    
[29] gplots_3.2.0                      ggnomics_0.1.2                   
[31] ggforce_0.3.2                     ggplot2_3.3.2                    
[33] gdata_2.18.0                      fuzzyjoin_0.1.6                  
[35] exact2x2_1.7.0                    exactci_1.4-5                    
[37] testthat_2.3.2                    ssanv_1.1                        
[39] dplyr_1.0.1                       doMC_1.3.6                       
[41] iterators_1.0.12                  foreach_1.5.0                    
[43] data.table_1.12.8                 copynumber_1.28.0                
[45] ComplexHeatmap_2.4.2              BSgenome.Hsapiens.UCSC.hg19_1.4.3
[47] BSgenome_1.56.0                   rtracklayer_1.48.0               
[49] Biostrings_2.56.0                 XVector_0.28.0                   
[51] GenomicRanges_1.40.0              GenomeInfoDb_1.24.0              
[53] IRanges_2.22.1                    S4Vectors_0.26.0                 
[55] BiocGenerics_0.34.0               rmarkdown_2.3                    

loaded via a namespace (and not attached):
  [1] tidyselect_1.1.0            lme4_1.1-23                
  [3] RSQLite_2.2.0               AnnotationDbi_1.68.0       
  [5] htmlwidgets_1.5.1           BiocParallel_1.22.0        
  [7] munsell_0.5.0               codetools_0.2-16           
  [9] statmod_1.4.34              withr_2.2.0                
 [11] colorspace_1.4-1            Biobase_2.48.0             
 [13] highr_0.8                   rstudioapi_0.11            
 [15] labeling_0.3                GenomeInfoDbData_1.2.3     
 [17] KMsurv_0.1-5                polyclip_1.10-0            
 [19] bit64_4.0.2                 farver_2.0.3               
 [21] vctrs_0.3.2                 generics_0.0.2             
 [23] xfun_0.16                   biovizBase_1.36.0          
 [25] BiocFileCache_1.12.0        R6_2.4.1                   
 [27] clue_0.3-57                 AnnotationFilter_1.12.0    
 [29] bitops_1.0-6                cachem_1.0.6               
 [31] DelayedArray_0.14.0         assertthat_0.2.1           
 [33] scales_1.1.1                nnet_7.3-14                
 [35] gtable_0.3.0                formula.tools_1.7.1        
 [37] ensembldb_2.12.1            rlang_1.0.6                
 [39] GlobalOptions_0.1.2         splines_4.0.2              
 [41] lazyeval_0.2.2              acepack_1.4.1              
 [43] dichromat_2.0-0             broom_0.7.0                
 [45] checkmate_2.0.0             modelr_0.1.8               
 [47] yaml_2.2.1                  GenomicFeatures_1.40.0     
 [49] backports_1.1.8             Hmisc_4.4-0                
 [51] tools_4.0.2                 ellipsis_0.3.1             
 [53] DNAcopy_1.64.0              plyr_1.8.6                 
 [55] Rcpp_1.0.7                  base64enc_0.1-3            
 [57] progress_1.2.2              zlibbioc_1.34.0            
 [59] RCurl_1.98-1.2              prettyunits_1.1.1          
 [61] rpart_4.1-15                openssl_1.4.2              
 [63] GetoptLong_1.0.2            cowplot_1.1.1              
 [65] zoo_1.8-8                   haven_2.3.1                
 [67] SummarizedExperiment_1.18.1 cluster_2.1.0              
 [69] fs_1.5.0                    circlize_0.4.10            
 [71] reprex_0.3.0                mitml_0.4-5                
 [73] ProtGenerics_1.20.0         matrixStats_0.56.0         
 [75] xtable_1.8-4                hms_0.5.3                  
 [77] evaluate_0.14               XML_3.99-0.3               
 [79] jpeg_0.1-8.1                readxl_1.3.1               
 [81] shape_1.4.4                 compiler_4.0.2             
 [83] biomaRt_2.44.0              mice_3.19.0                
 [85] KernSmooth_2.23-18          crayon_1.3.4               
 [87] minqa_1.2.4                 mgcv_1.8-31                
 [89] Formula_1.2-3               DBI_1.1.0                  
 [91] tweenr_1.0.1                dbplyr_1.4.4               
 [93] MASS_7.3-51.6               rappdirs_0.3.1             
 [95] boot_1.3-25                 Matrix_1.2-18              
 [97] cli_3.4.1                   pan_1.9                    
 [99] km.ci_0.5-6                 pkgconfig_2.0.3            
[101] GenomicAlignments_1.24.0    foreign_0.8-81             
[103] xml2_1.3.2                  rvest_0.3.6                
[105] VariantAnnotation_1.34.0    digest_0.6.25              
[107] cellranger_1.1.0            survMisc_0.5.6             
[109] htmlTable_2.0.1             operator.tools_1.6.3       
[111] curl_4.3                    Rsamtools_2.4.0            
[113] gtools_3.8.2                jomo_2.7-6                 
[115] rjson_0.2.20                nloptr_1.2.2.2             
[117] lifecycle_0.2.0             nlme_3.1-148               
[119] askpass_1.1                 pillar_1.4.6               
[121] lattice_0.20-41             KEGGREST_1.46.0            
[123] fastmap_1.0.1               httr_1.4.2                 
[125] glue_1.4.1                  png_0.1-7                  
[127] glmnet_4.1-10               bit_4.0.3                  
[129] stringi_1.4.6               blob_1.2.1                 
[131] latticeExtra_0.6-29         caTools_1.18.3             
[133] memoise_2.0.1
```

