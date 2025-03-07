# README for scRNA-seq and snRNA-seq Data Analysis Pipeline
This repository contains the R code for analyzing matrisome gene expression patterns in human brain development using single-cell RNA sequencing.

Preprint: Gim *et al.*, 2025 *"Deciphering Cell-Type and Temporal-Specific Matrisome Expression Signatures in Human Cortical Development and Neurodevelopmental Disorders via scRNA-Seq Meta-Analysis"*

https://www.biorxiv.org/content/10.1101/2025.02.24.639826v1.full

Extracellular matrix (ECM) plays a crucial role in guiding these processes, yet its specific contributions and the implications of its dysregulation in neurodevelopmental disorders (NDDs) remain underexplored.

![alt text](https://github.com/KANGBERGLAB/matrisome-meta-analysis/blob/dhg2022-patch-2/cortical_dev_ECM.png?raw=true)

The matrisome refers to the set of genes and proteins that compose and regulate the ECM. We found that 17.2% of core matrisome genes and 9.8% of matrisome-associated genes are reported as NDD risk genes

![alt text](https://github.com/KANGBERGLAB/matrisome-meta-analysis/blob/dhg2022-patch-2/fig2.png?raw=true)

In this study, we conducted a meta-analysis of single-cell RNA sequencing (scRNA-seq) data from 37 donors, gestational weeks (GWs) 8 to 26 across six independent studies to elucidate cell type-specific matrisome gene expression signatures and their dynamics in the developing human cortex.

![alt text](https://github.com/KANGBERGLAB/matrisome-meta-analysis/blob/dhg2022-patch-2/fig1.png?raw=true)


# System requirements

**Software dependencies**

R (≥ 4.1.0)

RStudio (≥ 2022.07.1)

Required R packages:
General packages: BiocManager, devtools, purrr, tidyverse, dplyr, tidyr, stringr, reshape2, ggplot2, cowplot, scales, ggrepel, ggpubr, ggprism, rstatix, ggsignif, plotly, ggthemes, RColorBrewer, wesanderson

Single-cell analysis: Seurat (≥ 4.3.0), SeuratWrappers, SeuratDisk, scuttle, DropletUtils, scran, uwot

Differential expression: DESeq2 (≥ 1.34.0), MAST, multtest, metap

Visualization: patchwork, ComplexHeatmap, colorRamp2, ggraph, VennDiagram, networkD3, gg3D, plot3D

Genomics: biomaRt, annotables, org.Hs.eg.db, EnsDb.Hsapiens.v86, HGNChelper

Network analysis: igraph, WGCNA, IReNA, monocle3, GENIE3

Regulatory analysis: BSgenome.Hsapiens.UCSC.hg38, TFBSTools, motifmatchr, RcisTarget


## Operating systems

macOS Monterey and later

Windows 10/11

**Tested configurations**

MacOS with M2 Max CPU and 96GB RAM

Windows 10/11 with vPro i9 and 128GB RAM

Hardware requirements
Minimum 32GB RAM (64GB+ recommended)

100GB free disk space

Multi-core processor (8+ cores recommended)



# Installation guide
**R and RStudio installation**
1. Download and install R from CRAN
https://cran.r-project.org/

2. Download and install RStudio from RStudio website
https://posit.co/download/rstudio-desktop/

3. Download matrisome DB
*Naba A, Clauser KR, Hoersch S, Liu H, Carr SA, and Hynes RO. The matrisome: in-silico definition and in-vivo characterization by proteomics of normal and tumor extracellular matrices. Molecular and Cellular Proteomics, 2012, 11(4):M111.014647. https://doi.org/10.1074/mcp.M111.014647*
https://sites.google.com/uic.edu/matrisome/matrisome-annotations/homo-sapiens

4. Custom functions
```
KANGBERGLAB/matrisome-meta-analysis
├── R/
│   ├── perform_deseq2/
│   ├── calculate_LISI_scores/
│   ├── create_volcano_plot/
│   ├── analyze_cell_type_degs/
│   └── analyze_temporal_patterns/
└── scripts/
    ├── ....Rmd
    └── ....Rmd
```


**Package installation**
```{r}
# Install Bioconductor and required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Seurat
install.packages("Seurat")

# Install other CRAN packages
install.packages(c("devtools", "tidyverse", "cowplot", "ggplot2", "patchwork", "igraph", 
                  "reshape2", "dplyr", "stringr", "tidyr", "ggpubr", "ggprism", 
                  "rstatix", "ggsignif", "HGNChelper", "scales", "VennDiagram", 
                  "plotly", "networkD3", "ggrepel", "RColorBrewer", "wesanderson", 
                  "ggthemes", "gg3D", "plot3D", "purrr", "WGCNA"))

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "MAST", "scuttle", "DropletUtils", "scran", 
                      "multtest", "metap", "biomaRt", "annotables", "org.Hs.eg.db", 
                      "EnsDb.Hsapiens.v86", "ComplexHeatmap", "BSgenome.Hsapiens.UCSC.hg38", 
                      "TFBSTools", "motifmatchr"))

# Install GitHub packages
devtools::install_github("satijalab/seurat-wrappers")
devtools::install_github("satijalab/seurat-disk")
devtools::install_github("aertslab/GENIE3")
devtools::install_github("cole-trapnell-lab/monocle3")
devtools::install_github("aertslab/RcisTarget")
```

# Demo
Expected run time for demo
30-40 minutes on a standard desktop computer

## Running the Demo and iLISI
Expected output
UMAP visualization showing cell type clustering, a data frame containing LISI, and mean iLISI and cLISI
```{r}
# Load required libraries
lapply(c("BiocManager","devtools","Seurat", "Matrix", "patchwork", "ggplot2", "ggthemes","scuttle", "DropletUtils", "uwot", "ggraph", "reshape2", "scran", "dplyr", "stringr", "ComplexHeatmap", "tidyr", "multtest", "metap", "generics", "igraph", 'ggpubr', 'ggprism', 'rstatix',"biomaRt","annotables","org.Hs.eg.db", "EnsDb.Hsapiens.v86", "tidyverse", 'ggsignif', 'HGNChelper', 'cowplot', 'colorRamp2', 'scales',"VennDiagram","plotly", "networkD3", "ggrepel", "RColorBrewer", "wesanderson", "DESeq2", "MAST", "SeuratDisk", "SeuratWrappers","gg3D","plot3D","purrr", "GO.db", "WGCNA","slingshot", "monocle3", "IReNA", "GENIE3","lisi","motifmatchr"), library, character.only = T)

# Load pre-processed Seurat object
load("demo_seurat.RData")

# Visualize the UMAP projection colored by cell type
DimPlot(merged_seurat, reduction = "umap", group.by = "sctype", 
        cols = c("#C7EBB1", "#FFD5FF", "#34854D", "#DEB0B0","#F26378", 
                 "#CA2F40", "#81B1DE","#C58862", "#7F517D", "#33608C", 
                 "#F5CA61", "#49525E"),
        label = TRUE, repel = TRUE) + NoAxes() + ggtitle("")

# Calculate iLISI and cLISI
results <- calculate_LISI_scores(
  object = merged_seurat,
  batch = "study",
  cluster.anno = "sctype",
  sample.size = 20000
)

cell_lisi_scores <- lisi_results$lisi_scores
mean_iLISI <- lisi_results$mean_iLISI
mean_cLISI <- lisi_results$mean_cLISI
print(mean_iLISI)
print(mean_cLISI)
```

## Performing Pseudobulk Differential Gene Expression (DGE) analysis
Expected output
A table of differentially expressed genes with log2fold change and adjusted p-value for each cell type
```{r}
# Duplicate object to avoid dataset contamination
merged_seurat2 <- merged_seurat
merged_seurat2[["RNA"]] <- JoinLayers(merged_seurat2[["RNA"]])
pseudo_obj <- AggregateExpression(merged_seurat, assays = "RNA", return.seurat = T, group.by = c("sctype", "donor"))
table(pseudo_obj$sctype)
Idents(pseudo_obj) <- pseudo_obj$sctype

colpal <- c ("#C7EBB1", "#FFD5FF", "#34854D", "#DEB0B0","#F26378", "#CA2F40", "#81B1DE","#C58862", "#7F517D", "#33608C", "#F5CA61", "#49525E")
cell_levels <- c("vRG","oRG","tRG", "IPC", "Imm.GExN","M.GExN","Unknown","GABA IN","Astrocyte","OPC","Microglia","Endothelial")

pseudo_obj@active.ident <- factor(x = pseudo_obj@active.ident, levels = cell_levels)
pseudo_obj <- NormalizeData(pseudo_obj)
pseudo_obj <- ScaleData(pseudo_obj,features = rownames(pseudo_obj))

# Check donors per cell type
n_per_celltype <- table(pseudo_obj$sctype)
n_per_celltype



# Get all unique cell types
cell_types <- unique(Idents(pseudo_obj))

# Perform DESeq2 for each cell type
all_results <- lapply(cell_types, function(ct) {
  result <- perform_deseq2(pseudo_obj, ct)
  result$cell_type <- ct
  return(result)
})

# Combine all results
final_results <- do.call(rbind, all_results)
head(final_results)

# Save results
write.csv(final_results, file = "~/your/directory/all_cell_types_de_results.csv", row.names = FALSE)
```

## Temporal dynamics
Expected output
Top temporally increasing and decreasing genes (ranked by +/-gradient) with Pearson's R and p-values for each cell type. Visualization of linear plots per cell type.
```{r}
# Set project directory
project_dir <-  "~/your/directory"

# Run analysis
results <- analyze_temporal_patterns (
  merged_seurat2,
  project_dir,
  split.by = c("NSC", "IPC", "Imm.GExN", "M.GExN",
                      "Unknown", "GABA IN", "Astrocyte", "OPC", 
                      "Microglia", "Endothelial"),
  trend.in = "trimester",
  trend.order = c("first", "early", "mid"),
  pseudobulk.by = "donor",
  suppress.warnings = TRUE
)

# View results
print(results$plot)
print(results$rankings)
```
# Instructions for use
## Data preparation
1. The full analysis requires integration of multiple scRNA-seq datasets from different studies of brain development.
2. Download raw data from GEO (accession numbers provided in the publication).
3. Download functions and scripts from this repository.
4. Save the raw data in a directory structure as follows:
```
project_dir/
├── raw_data/
│   ├── bhaduri/
│   ├── bruggen/
│   ├── cameron/
│   ├── eze/
│   ├── trevino/
│   └── zhong/
└── scripts/
```
## Performing analyses
**data_integration.Rmd**: Data integration and preprocessing

This script imports individual datasets, performs quality control, normalizes, and integrates data.

Due to the large memory requirements, it's recommended to run this on a high-memory machine.

Run time: 4-6 hours with 96GB RAM

**other analysis**: Matrisome-focused analysis

This script performs comprehensive analysis of matrisome gene expression in the integrated dataset.

It contains 17 separate analysis sections that can be run independently.

Run time for complete analysis: 8-12 hours

**Output interpretation**
The analyses generate various output files including:

• Annotated Seurat object

• Differential expression results for matrisome genes

• Cell type-specific matrisome markers

• Regulatory network files

• High-quality visualizations

**Memory considerations**
The analysis was designed and tested on high-memory configurations. If running on systems with less than 64GB RAM:

Consider using a subset of cells (random sampling)

Run sections individually rather than the entire Rmd files

Use functions like gc() frequently to free memory

Process and save intermediate results to disk

**Customization**
The scripts can be modified to focus on specific cell types or developmental stages by adjusting the cell subsetting parameters. The matrisome gene list can also be replaced with custom gene sets to study other pathways or processes.

## Summary of ey analyses:
1. Dataset loading, QC, normalization, and batch effect correction
2. Cell clustering and UMAP visualization
3. Integration quality assessment
4. Cell type annotation using reference markers
5. Differential expression analysis with pseudobulk approach
6. Identification of cell type-specific matrisome signatures
7. Characterization of neural stem cell subtypes
8. Co-expression and regulatory network analysis
9. Temporal expression analysis
10. Cell-cell interaction inference

# Session information
sessionInfo()
```
R version 4.4.2 (2024-10-31)
Platform: x86_64-apple-darwin20
Running under: macOS Sequoia 15.3.1

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] future_1.34.0               CellChat_2.1.2              GENIE3_1.26.0               IReNA_0.99.0               
 [5] RcisTarget_1.24.0           DDRTree_0.1.5               irlba_2.3.5.1               pheatmap_1.0.12            
 [9] monocle3_1.3.7              slingshot_2.12.0            TrajectoryUtils_1.12.0      princurve_2.1.6            
[13] WGCNA_1.73                  fastcluster_1.2.6           dynamicTreeCut_1.63-1       GO.db_3.19.1               
[17] plot3D_1.4.1                gg3D_0.0.0.9000             SeuratWrappers_0.3.5        SeuratDisk_0.0.0.9021      
[21] MAST_1.30.0                 DESeq2_1.44.0               wesanderson_0.3.7           RColorBrewer_1.1-3         
[25] ggrepel_0.9.6               networkD3_0.4               plotly_4.10.4               VennDiagram_1.7.3          
[29] futile.logger_1.4.3         scales_1.3.0                colorRamp2_0.1.0            cowplot_1.1.3              
[33] HGNChelper_0.8.15           ggsignif_0.6.4              lubridate_1.9.3             forcats_1.0.0              
[37] purrr_1.0.2                 readr_2.1.5                 tibble_3.2.1                tidyverse_2.0.0            
[41] EnsDb.Hsapiens.v86_2.99.0   ensembldb_2.28.1            AnnotationFilter_1.28.0     GenomicFeatures_1.56.0     
[45] org.Hs.eg.db_3.19.1         AnnotationDbi_1.66.0        annotables_0.2.0            biomaRt_2.60.1             
[49] rstatix_0.7.2               ggprism_1.0.5               ggpubr_0.6.0                igraph_2.1.2               
[53] generics_0.1.3              metap_1.11                  multtest_2.60.0             tidyr_1.3.1                
[57] ComplexHeatmap_2.20.0       stringr_1.5.1               dplyr_1.1.4                 scran_1.32.0               
[61] reshape2_1.4.4              ggraph_2.2.1                uwot_0.2.2                  DropletUtils_1.24.0        
[65] scuttle_1.14.0              SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0 Biobase_2.64.0             
[69] GenomicRanges_1.56.2        GenomeInfoDb_1.40.1         IRanges_2.38.1              S4Vectors_0.42.1           
[73] BiocGenerics_0.50.0         MatrixGenerics_1.16.0       matrixStats_1.4.1           ggthemes_5.1.0             
[77] ggplot2_3.5.1               patchwork_1.3.0             Matrix_1.7-1                Seurat_5.1.0               
[81] devtools_2.4.5              usethis_3.1.0               BiocManager_1.30.25         SeuratObject_5.0.2         
[85] sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] IRkernel_1.3.2                          DBI_1.2.3                               bslib_0.8.0                            
  [4] httr_1.4.7                              registry_0.5-1                          BiocParallel_1.38.0                    
  [7] prettyunits_1.2.0                       yulab.utils_0.1.8                       ProtGenerics_1.36.0                    
 [10] ggplotify_0.1.2                         GenomicAlignments_1.40.0                sparseMatrixStats_1.16.0               
 [13] spatstat.geom_3.3-4                     splitstackshape_1.4.8                   pillar_1.9.0                           
 [16] R6_2.5.1                                boot_1.3-31                             mime_0.12                              
 [19] reticulate_1.40.0                       edgeR_4.2.2                             monocle_2.32.0                         
 [22] viridis_0.6.5                           Rhdf5lib_1.26.0                         ROCR_1.0-11                            
 [25] Hmisc_5.2-1                             rprojroot_2.0.4                         limma_3.60.6                           
 [28] parallelly_1.40.1                       GlobalOptions_0.1.2                     FNN_1.1.4.1                            
 [31] rbibutils_2.3                           caTools_1.18.3                          polyclip_1.10-7                        
 [34] NMF_0.28                                beachmat_2.20.0                         htmltools_0.5.8.1                      
 [37] fansi_1.0.6                             lambda.r_1.2.4                          remotes_2.5.0                          
 [40] fastICA_1.2-5.1                         car_3.1-3                               ChIPseeker_1.40.0                      
 [43] fgsea_1.30.0                            spatstat.utils_3.1-1                    clusterProfiler_4.12.6                 
 [46] TFisher_0.2.0                           rpart_4.1.23                            clue_0.3-66                            
 [49] scatterpie_0.2.4                        fitdistrplus_1.2-1                      goftest_1.2-3                          
 [52] tidyselect_1.2.1                        RSQLite_2.3.9                           GenomeInfoDbData_1.2.12                
 [55] utf8_1.2.4                              ScaledMatrix_1.12.0                     scattermore_1.2                        
 [58] sessioninfo_1.2.2                       spatstat.data_3.1-4                     gridExtra_2.3                          
 [61] fs_1.6.5                                sctransform_0.4.1                       future.apply_1.11.3                    
 [64] graph_1.82.0                            R.oo_1.27.0                             RcppHNSW_0.6.0                         
 [67] mathjaxr_1.6-0                          uuid_1.2-1                              rtracklayer_1.64.0                     
 [70] Rtsne_0.17                              DelayedMatrixStats_1.26.0               lazyeval_0.2.2                         
 [73] sass_0.4.9                              carData_3.0-5                           munsell_0.5.1                          
 [76] treeio_1.28.0                           R.utils_2.12.3                          profvis_0.4.0                          
 [79] bitops_1.0-9                            R.methodsS3_1.8.2                       labeling_0.4.3                         
 [82] KEGGREST_1.44.1                         promises_1.3.2                          shape_1.4.6.1                          
 [85] rhdf5filters_1.16.0                     zoo_1.8-12                              locfit_1.5-9.10                        
 [88] DelayedArray_0.30.1                     arrow_18.1.0                            RSpectra_0.16-2                        
 [91] multcomp_1.4-26                         assertthat_0.2.1                        tools_4.4.2                            
 [94] ape_5.8                                 shiny_1.9.1                             BiocFileCache_2.12.0                   
 [97] rlang_1.1.4                             BiocSingular_1.20.0                     ggridges_0.5.6                         
[100] evaluate_1.0.1                          httr2_1.0.7                             BiocIO_1.14.0                          
[103] plotrix_3.8-4                           colorspace_2.1-1                        ellipsis_0.3.2                         
[106] data.table_1.16.4                       withr_3.0.2                             presto_1.0.0                           
[109] RCurl_1.98-1.16                         restfulr_0.0.15                         xtable_1.8-4                           
[112] plyr_1.8.9                              lme4_1.1-35.5                           aplot_0.2.3                            
[115] systemfonts_1.1.0                       httpuv_1.6.15                           rmarkdown_2.29                         
[118] metapod_1.12.0                          RCy3_2.24.0                             MASS_7.3-61                            
[121] dqrng_0.4.1                             broom_1.0.7                             deldir_2.0-4                           
[124] sandwich_3.1-1                          rhdf5_2.48.0                            tensor_1.5                             
[127] vctrs_0.6.5                             lifecycle_1.0.4                         codetools_0.2-20                       
[130] fastDummies_1.7.4                       mnormt_2.1.1                            here_1.0.1                             
[133] nlme_3.1-166                            combinat_0.0-8                          progress_1.2.3                         
[136] dbplyr_2.5.0                            jquerylib_0.1.4                         pkgload_1.4.0                          
[139] Rcpp_1.0.13-1                           rstudioapi_0.17.1                       stringi_1.8.4                          
[142] VGAM_1.1-12                             hms_1.1.3                               pbapply_1.7-2                          
[145] minqa_1.2.8                             cachem_1.1.0                            base64url_1.4                          
[148] tcltk_4.4.2                             hdf5r_1.3.11                            tidytree_0.4.6                         
[151] listenv_0.9.1                           XVector_0.44.0                          urlchecker_1.0.1                       
[154] mutoss_0.1-13                           ggtree_3.12.0                           enrichplot_1.24.4                      
[157] GetoptLong_1.0.5                        pkgbuild_1.4.5                          ggfun_0.1.8                            
[160] HDF5Array_1.32.1                        IRdisplay_1.1                           misc3d_0.9-1                           
[163] SparseArray_1.4.8                       htmlwidgets_1.6.4                       Formula_1.2-5                          
[166] memoise_2.0.1                           crayon_1.5.3                            gridGraphics_0.5-1                     
[169] rappdirs_0.3.3                          S4Arrays_1.4.1                          xml2_1.3.6                             
[172] filelock_1.0.3                          preprocessCore_1.66.0                   GOSemSim_2.30.2                        
[175] UCSC.utils_1.0.0                        png_0.1-8                               progressr_0.15.1                       
[178] tzdb_0.4.0                              fastmap_1.2.0                           GSEABase_1.66.0                        
[181] coda_0.19-4.1                           tidygraph_1.3.1                         pkgconfig_2.0.3                        
[184] cli_3.6.3                               DOSE_3.30.5                             ggforce_0.4.2                          
[187] nnet_7.3-19                             gridBase_0.4-7                          ggalluvial_0.12.5                      
[190] lmtest_0.9-40                           RcppAnnoy_0.0.22                        slam_0.1-55                            
[193] timechange_0.3.0                        viridisLite_0.4.2                       foreign_0.8-87                         
[196] splines_4.4.2                           blob_1.2.4                              impute_1.78.0                          
[199] annotate_1.82.0                         XML_3.99-0.17                           network_1.18.2                         
[202] numDeriv_2016.8-1.1                     globals_0.16.3                          knitr_1.49                             
[205] ica_1.0-3                               spam_2.11-0                             compiler_4.4.2                         
[208] rjson_0.2.23                            qqconf_1.3.2                            bit_4.5.0.1                            
[211] leidenbase_0.1.31                       sn_2.1.1                                BiocNeighbors_1.22.0                   
[214] glue_1.8.0                              ggnetwork_0.5.13                        formatR_1.14                           
[217] digest_0.6.37                           leiden_0.4.3.1                          pbdZMQ_0.3-13                          
[220] graphlayouts_1.2.1                      foreach_1.5.2                           spatstat.random_3.3-2                  
[223] zlibbioc_1.50.0                         dotCall64_1.2                           tweenr_2.0.3                           
[226] lattice_0.22-6                          statmod_1.5.0                           rsvd_1.0.5                             
[229] gson_0.1.0                              nloptr_2.1.1                            mvtnorm_1.3-2                          
[232] yaml_2.3.10                             qvalue_2.36.0                           later_1.4.1                            
[235] statnet.common_4.10.0                   backports_1.5.0                         shadowtext_0.1.4                       
[238] Rsamtools_2.20.0                        parallel_4.4.2                          sna_2.8                                
[241] miniUI_0.1.1.1                          gtable_0.3.6                            abind_1.4-8                            
[244] xfun_0.49                               Biostrings_2.72.1                       curl_6.0.1                             
[247] doParallel_1.0.17                       KernSmooth_2.23-24                      futile.options_1.0.1                   
[250] survival_3.7-0                          jsonlite_1.8.9                          magrittr_2.0.3                         
[253] svglite_2.1.3                           base64enc_0.1-3                         iterators_1.0.14                       
[256] TH.data_1.1-2                           spatstat.univar_3.1-1                   HSMMSingleCell_1.24.0                  
[259] fastmatch_1.1-4                         checkmate_2.3.2                         gtools_3.9.5                           
[262] htmlTable_2.4.3                         spatstat.sparse_3.1-0                   AUCell_1.26.0                          
[265] rngtools_1.5.2                          RANN_2.6.2                              bluster_1.14.0                         
[268] repr_1.1.7                              circlize_0.4.16                         TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
[271] spatstat.explore_3.3-3                  RJSONIO_1.3-1.9                         bit64_4.5.2                            
[274] cluster_2.1.7                           Rdpack_2.6.2                            farver_2.1.2                           
[277] gplots_3.2.0        
```
