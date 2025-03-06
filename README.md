README for scRNA-seq and snRNA-seq Data Analysis Pipeline
1. System Requirements
Software Dependencies
•	Operating System: Linux (Ubuntu 20.04+), macOS (Monterey+), Windows 10+
•	R Version: 4.4.2
•	RStudio Version: 2023.06+
•	Required R Packages: 
o	Seurat v5.1.0
o	tidyverse v2.0.0
o	dplyr v1.1.4
o	DESeq2 v1.44.0
o	ComplexHeatmap v2.20.0
o	ggplot2 v3.5.1
o	WGCNA v1.73
o	monocle3 v1.3.7
o	LISI v1.0.0
o	ggpubr v0.6.0
o	CellChat v1.6.1
o	motifmatchr v1.28.0
•	Additional Requirements: 
o	Python (for some auxiliary scripts): Python 3.8+
o	Hardware: No special hardware required, but a system with 16GB+ RAM is recommended for large datasets.
Versions Tested
•	Windows 10
•	R v4.4.2
•	RStudio 2023.06
2. Installation Guide
Installation Steps
1.	Install R and RStudio from CRAN and RStudio.
2.	Install required R packages using the following command in R: 
3.	install.packages(c("Seurat", "tidyverse", "dplyr", "DESeq2", "ComplexHeatmap", "ggplot2", "WGCNA", "monocle3", "LISI", "ggpubr", "CellChat", "motifmatchr"))
4.	(Optional) Install Python dependencies if using additional scripts: 
5.	pip install numpy pandas
Typical Installation Time
•	On a standard desktop (16GB RAM, Intel i7 processor): ~15 minutes
3. Demo
Running the Demo
1.	Open RStudio and set the working directory to the project folder.
2.	Load the necessary libraries: 
3.	library(Seurat)
4.	library(tidyverse)
5.	library(dplyr)
6.	Load the demo dataset: 
7.	demo_data <- readRDS("demo_dataset.rds")
8.	Process and integrate data: 
9.	demo_data <- NormalizeData(demo_data)
10.	demo_data <- FindVariableFeatures(demo_data)
11.	demo_data <- ScaleData(demo_data)
12.	Perform clustering and visualization: 
13.	demo_data <- RunPCA(demo_data)
14.	demo_data <- FindNeighbors(demo_data, dims = 1:20)
15.	demo_data <- FindClusters(demo_data, resolution = 1.5)
16.	demo_data <- RunUMAP(demo_data, dims = 1:20)
17.	DimPlot(demo_data, reduction = "umap")
Expected Output
•	A UMAP plot of clustered single-cell RNA-seq data.
•	Identification of distinct clusters with color-coded visualization.
Expected Run Time
•	On a standard desktop (16GB RAM, Intel i7 processor): ~20 minutes
4. Instructions for Use
Running the Software on Custom Data
1.	Load the count matrix: 
2.	count_matrix <- Read10X(data.dir = "path_to_your_data")
3.	Create a Seurat object: 
4.	seurat_obj <- CreateSeuratObject(counts = count_matrix, min.cells = 3, min.features = 200)
5.	Perform quality control: 
6.	seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
7.	Normalize data and identify variable features: 
8.	seurat_obj <- NormalizeData(seurat_obj)
9.	seurat_obj <- FindVariableFeatures(seurat_obj)
10.	Run PCA and clustering: 
11.	seurat_obj <- ScaleData(seurat_obj)
12.	seurat_obj <- RunPCA(seurat_obj)
13.	seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
14.	seurat_obj <- FindClusters(seurat_obj, resolution = 1.5)
15.	Visualize using UMAP: 
16.	seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
17.	DimPlot(seurat_obj, reduction = "umap")
Additional Functionalities
•	Differential Gene Expression Analysis: 
•	de_results <- FindMarkers(seurat_obj, ident.1 = "Cluster1", ident.2 = "Cluster2", min.pct = 0.25)
•	Cell Type Annotation using scType: 
•	cell_types <- scType(seurat_obj, marker_gene_table = "Supplementary_Table_1.csv")
•	Integration of Multiple Datasets: 
•	integrated_data <- IntegrateData(anchorset = integration_anchors)
Output Files
•	Clustered scRNA-seq data with annotations (RDS file)
•	UMAP plot (PDF/PNG file)
•	Differential expression results (CSV file)
