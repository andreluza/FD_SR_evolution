Could evolution explain differences in functional richness and evenness
along species richness gradients? A formal test using multivariate trait
simulations
================
2026-02-20

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

#### Organization of the repository: 

*Root*  
\|  
\|– *data_fish*: data of reef fishes from Morais et al. 2017 - file
“UpdatedData_RMorais_et_al_2017.csv”.  \|——– The folder ‘TACT’ hosts the
fish phylogenies (from Siqueira et al. 2020), file
“Reef_fish_all_combined.trees”.  
\|——– The file
“Atributos_especies_Atlantico\_&\_Pacifico_Oriental_2020_04_28.csv” has
trait data from Quimbayo et al. 2021.  
\|– *data rodents*:  
\|——– “AppendixS1- Small_mammal_data.csv” = data from Luza et al. 2019.
“Mammal_Communities.csv” = data from Figueiredo et al. 2017
(“Localities.csv” and “References.csv” are metadata associated to data
of Figueiredo et al. 2017).  
\|——– “limites_integradores_wgs84_v1_2_0” folder with Atlantic forest
shapefiles Muylaert et al. 2018.  
\|——– “Penone_et_al_2016_mammal_trait_data_imputed.csv” = mammal trait
data from Penone et al. 2016.  
\|——– “Sigmodontinae_413species100Trees.trees” = rodent phylogenies from
Upham et al. 2016.  
\|– *Processed data*  
\|——– Basic data to fit models (“image*.RData”) for each group. Data of
environmental and historical covariates after extraction to survey sites
are available in this folder.  
\|– *Output*: files in format ‘RData’ with parameters of
macroevolutionary models fitted to the data (sufix ’simulated*.RData’)
and Bayesian linear models fitted to FRic and FEve data (sufix”GLM”).  
\|——– Figures: folder with figures shown in the main text and supporting
information.  
\|– *R*: R scripts used in the analyzes. The codes are numbered from 1
to 15, depicting different steps of analysis.  

<!-- badges: start -->
<!-- badges: end -->

#### This paper was produced using the following software and associated packages: 

    ## R version 4.5.2 (2025-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=Portuguese_Brazil.utf8  LC_CTYPE=Portuguese_Brazil.utf8   
    ## [3] LC_MONETARY=Portuguese_Brazil.utf8 LC_NUMERIC=C                      
    ## [5] LC_TIME=Portuguese_Brazil.utf8    
    ## 
    ## time zone: Europe/Paris
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] worrms_0.4.3            terra_1.8-86            sf_1.0-23              
    ##  [4] scatterpie_0.2.6        ggrepel_0.9.6           rnaturalearthdata_1.0.0
    ##  [7] rnaturalearth_1.1.0     sjlabelled_1.2.0        sjmisc_2.8.11          
    ## [10] sjPlot_2.9.0            cmdstanr_0.9.0          emmeans_2.0.1          
    ## [13] brms_2.23.0             Rcpp_1.1.0              tidybayes_3.0.7        
    ## [16] scales_1.4.0            gridExtra_2.3           ggplot2_4.0.1          
    ## [19] fishtree_0.3.4          FD_1.0-12.3             geometry_0.5.2         
    ## [22] ade4_1.7-23             missForest_1.6.1        mvMORPH_1.2.1          
    ## [25] subplex_1.9             corpcor_1.6.10          picante_1.8.2          
    ## [28] nlme_3.1-168            vegan_2.7-2             permute_0.9-8          
    ## [31] geiger_2.0.11           phytools_2.5-2          maps_3.4.3             
    ## [34] ape_5.8-1               reshape_0.8.10          here_1.0.2             
    ## [37] openxlsx_4.2.8.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      tensorA_0.36.2.1        rstudioapi_0.17.1      
    ##   [4] jsonlite_2.0.0          magrittr_2.0.4          estimability_1.5.1     
    ##   [7] farver_2.1.2            rmarkdown_2.30          fs_1.6.6               
    ##  [10] vctrs_0.6.5             memoise_2.0.1           htmltools_0.5.9        
    ##  [13] itertools_0.1-3         distributional_0.5.0    DEoptim_2.2-8          
    ##  [16] deSolve_1.40            KernSmooth_2.23-26      plyr_1.8.9             
    ##  [19] cachem_1.1.0            igraph_2.2.1            lifecycle_1.0.4        
    ##  [22] iterators_1.0.14        pkgconfig_2.0.3         Matrix_1.7-4           
    ##  [25] R6_2.6.1                fastmap_1.2.0           rbibutils_2.4          
    ##  [28] magic_1.6-1             digest_0.6.39           numDeriv_2016.8-1.1    
    ##  [31] ps_1.9.1                rprojroot_2.1.1         clusterGeneration_1.3.8
    ##  [34] randomForest_4.7-1.2    polyclip_1.10-7         abind_1.4-8            
    ##  [37] mgcv_1.9-3              compiler_4.5.2          proxy_0.4-29           
    ##  [40] rngtools_1.5.2          withr_3.0.2             doParallel_1.0.17      
    ##  [43] S7_0.2.1                backports_1.5.0         optimParallel_1.0-2    
    ##  [46] DBI_1.2.3               ggforce_0.5.0           MASS_7.3-65            
    ##  [49] rappdirs_0.3.3          classInt_0.4-11         scatterplot3d_0.3-44   
    ##  [52] loo_2.9.0               units_1.0-0             tools_4.5.2            
    ##  [55] ranger_0.17.0           otel_0.2.0              zip_2.3.3              
    ##  [58] glue_1.8.0              quadprog_1.5-8          grid_4.5.2             
    ##  [61] checkmate_2.3.3         cluster_2.1.8.1         generics_0.1.4         
    ##  [64] gtable_0.3.6            class_7.3-23            tidyr_1.3.2            
    ##  [67] foreach_1.5.2           pillar_1.11.1           ggdist_3.3.3           
    ##  [70] stringr_1.6.0           yulab.utils_0.2.3       spam_2.11-1            
    ##  [73] posterior_1.6.1         splines_4.5.2           tweenr_2.0.3           
    ##  [76] dplyr_1.1.4             lattice_0.22-7          tidyselect_1.2.1       
    ##  [79] knitr_1.51              arrayhelpers_1.1-0      xfun_0.55              
    ##  [82] expm_1.0-0              bridgesampling_1.2-1    matrixStats_1.5.0      
    ##  [85] stringi_1.8.7           ggfun_0.2.0             yaml_2.3.12            
    ##  [88] evaluate_1.0.5          codetools_0.2-20        tibble_3.3.0           
    ##  [91] cli_3.6.5               RcppParallel_5.1.11-1   xtable_1.8-4           
    ##  [94] pbmcapply_1.5.1         Rdpack_2.6.4            processx_3.8.6         
    ##  [97] glassoFast_1.0.1        coda_0.19-4.1           svUnit_1.0.8           
    ## [100] rstantools_2.5.0        dotCall64_1.2           doRNG_1.8.6.2          
    ## [103] bayesplot_1.15.0        Brobdingnag_1.2-9       phangorn_2.12.1        
    ## [106] mvtnorm_1.3-3           e1071_1.7-17            insight_1.4.4          
    ## [109] purrr_1.2.0             combinat_0.0-8          rlang_1.1.6            
    ## [112] fastmatch_1.1-6         mnormt_2.1.1
