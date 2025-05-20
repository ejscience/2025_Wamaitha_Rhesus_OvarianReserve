# Analysis of published cyno datasets -------------------------------
library(Seurat)
library(tidyverse)
library(DropletUtils)
library(BiocManager)
library(RSpectra)
library(hdf5r)
library(scDblFinder)
library(orthogene)
library(harmony)
library(patchwork)
library(viridis)
library(khroma)

-------------------------------
# Import Option 1 - From count matrix (.txt) - GSE149629-------------------------------
## Create Seurat object, identify each sample/replicate from barcode suffix data
cynoT2 <- read.table("~/Box/Genomics/10x datasets/Zhao_E84_E116/GSE149629_Count_matrix.txt.gz")

cyno <- CreateSeuratObject (counts = cynoT2, project = "Zhao2020", assay = "RNA", min.cells = 3, min.features = 200, names.field = 2)
### names.field in this instance is barcode_sampleno.
### Assign sample identities
table(Idents(cyno))
cyno <- RenameIdents(cyno, `1` = "Zhao_E84")
cyno <- RenameIdents(cyno, `2` = "Zhao_E116")
table(Idents(cyno))
cyno[["orig.ident"]] <- Idents(cyno)

saveRDS(cyno, "Zhao_T2.rds")


# Import Option 2 - From 10x export files (barcodes.tsv, features.tsv, matrix.mtx) - GSE194264, GSE160043 
--------------------------------
## Create individual Seurat objects from each sample/replicate
data <- Read10X(data.dir = "~/Box/Genomics/10x datasets/Mizuta_W18_filtered/")
W18 <- CreateSeuratObject(counts = data, project = "Mizuta_W18")
## Repeat for each sample/replicate

## Note: if there are errors, check for empty rows as follows - can replace with 0 for a quick workaround
any(rownames(data)=="")
which(rownames(data)=="")
rownames(data)[6114] <- 0


# Import Option 3 - From h5 files (.h5) - Zenodo 6918355  ------------------------------
## Create Seurat object, identify each sample/replicate from barcode suffix data
## For MacOS, install hdf5 dependencies in Terminal - first install Homebrew + Xcode command line tools, then MacPorts (https://www.macports.org/install.php)
sudo port install hdf5
## Install hdf5r reader, then load library
install.packages("hdf5r")

data <- Read10X_h5("~/Box/Genomics/10x datasets/Chen_E52/monkey_OV_E52_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
E52 <- CreateSeuratObject(counts = data, project = "Chen_E52")
## Repeat for each sample/replicate

## For options 2 and 3 - merge sample Seurat objects into one larger object dataset -------------------------------
cynoT2 <- merge(W8_1, y = c(W8_2, W10, W12, W16_1, W16_2, W18), add.cell.ids = c("W8_1","W8_2", "W10", "W12", "W16_1", "W16_2", "W18"), project = "Mizuta2022", merge.data = TRUE)

cynoT2
head(colnames(cynoT2))
tail(colnames(cynoT2))
unique(sapply(X = strsplit(colnames(cynoT2), split = "_"), FUN = "[", 1))
table(cynoT2$orig.ident)
cyno_merge[["sample"]] <- Idents(cyno_merge)

saveRDS(cynoT2, "Mizuta_T2.rds")
## Repeat for each sample/replicate to merge into datasets


# Merge all datasets ---------------------------------------------------
mT2 <- readRDS("Mizuta_T2.rds")
zT2 <- readRDS("Zhao_T2.rds")
cT1 <- readRDS("Chen_T1.rds")
sT1 <- readRDS("Sasaki_T1.rds")

cyno_combined <- merge(mT2, y = c(zT2, cT1, sT1), add.cell.ids = c("Mizuta", "Zhao", "Chen", "Sasaki"), project = "Cyno", merge.data = TRUE)

head(cyno_combined@meta.data
     
##Save file
saveRDS(cyno_combined, "Cyno_merge.rds")

------------------------------------     
# Pre-process data - QC mito + ribo, cell cycle ---------------------- 
## QC - percent mito, percent ribo
## Check mito genes
grep ("^COX3", rownames(cyno_combined[["RNA"]]),value = T)
### Note - "^MT-" / "^mt-" refers to human genome mitochondrial annotation; consequently, classic percent.mito command comes back as 0 across all samples - can use grep command to search for specific genes and add similar to percent.ribo command.
### Sasaki 2021 "Mitochondrial content was not used for normalization since the majority of mitochondrial genes starting with ‘‘MT-’’ were not annotated in MacFas5.0, a reference genome for cynomolgus monkeys." - there are no mito genes in these samples!
### Chen 2022 "In human the mitochondria-related genes were the genes starting with MT, and ND1,ND2, COX1, COX2, ATP8, ATP6, COX3, ND3, ND4L, ND4, ND5, ND6 and CYTB were the mitochondria-related genes used in other three species."
### Wamaitha cyno data combo: ND1, ND2, ND3, ND4L, ND4, ND5, ND6, COX3, and CYTB
cyno_combined[["percent.mt"]] <- PercentageFeatureSet(cyno_combined, pattern = "^ND[[:digit:]]|^CYTB|^COX3")
cyno_combined[["percent.ribo"]] <- PercentageFeatureSet(cyno_combined, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
     
plot(x = cyno_combined$nCount_RNA, y = cyno_combined$percent.mt, main = "UMI Counts x Percent Mito")
     
cyno_filter <- subset(cyno_combined, subset = percent.mt < 20)
     
## QC - cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cyno_filter <- CellCycleScoring(cyno_filter, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(cyno_filter)
     
table(Idents(cyno_filter))
Idents(cyno_filter) <- "orig.ident"
     
## QC - nFeature, nRNA filter
cyno_filter <- subset(cyno_filter, subset = nFeature_RNA > 200)
     
### Satija basic: pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
### Sasaki 2021 -Cells with more than 5500 genes or fewer than 1000 genes were filtered out 
### Zhao 2020 -Cells with more than 3000 genes or fewer than 200 genes were filtered out (in Garcia-Alonso 2022, cells with fewer than 300 genes)
     
## Check (and optional save) QC plots
Idents(cyno_filter) <- "orig.ident"
VlnPlot(cyno_filter, raster = FALSE, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.1, ncol = 2)

plot1 <- FeatureScatter(cyno_filter, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cyno_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
     
ggsave("cyno_filterStats.png", units="in", width = 8, height = 8, device='png', dpi=300)
     
## Save file 
saveRDS(cyno_filter, "Cyno_merge.rds")
 
--------------------------------------    
# QC - Exclude doublets 
## Identify doublets using scDblfinder -------------
## Tutorial at https://biostatsquid.com/scdblfinder-tutorial/
## Don't have to split by sample, unlike DoubletFinder or alternatives-----------------------------------
library(BiocManager)
BiocManager::install("scDblFinder")
library(scDblFinder)
project_path <- "E:/10X Data/Analysis_Sissy/10x datasets"
in_path <- project_path
     
## Load data
seu <- readRDS(file.path(in_path, "Cyno_merge.rds"))
head(seu@meta.data)
table(seu$orig.ident)
     
## Convert to single cell experiment + run scDblFinder
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce, samples = "orig.ident")
table(sce$scDblFinder.class)
sce@colData@listData %>% as.data.frame() %>% head()

## Add metadata to Seurat object
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>%
dplyr::select(starts_with('scDblFinder'))
head(meta_scdblfinder)
rownames(meta_scdblfinder) <- sce@colData@rownames
head(meta_scdblfinder)
seu <- AddMetaData(object = seu, metadata = meta_scdblfinder %>% dplyr::select('scDblFinder.class'))
table(seu$scDblFinder.class)
     
## Visualise singlet vs doublet metrics
VlnPlot(seu, raster = FALSE, group.by = 'orig.ident', split.by = "scDblFinder.class",
features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
 ncol =2, pt.size = 0)
     
## Exclude doublets
Idents(seu) <- "scDblFinder.class"
table(Idents(seu))
seu <- subset(x = seu, idents = 'singlet')
table(Idents(seu))
     
## Save file
saveRDS(seu, "Cyno_merge.rds")

----------------------------------------
# Convert cyno to rhesus orthologs: Orthogene
BiocManager::install("orthogene")
library(orthogene)

## Find species
species_dt <- all_species()
### gprofiler has both M.mulatta and M.fascicularis; homologene and babelgene only M.mulattta

cyno_merge <- readRDS(file.path(in_path, "Cyno_merge.rds"))

## Convert Seurat object to sparse matrix
data <- GetAssayData(cyno_merge[["RNA"]], layer = "data")

## Convert orthologs
gene_df <- orthogene::convert_orthologs(gene_df = data,
                                        gene_input = "rownames",
                                        gene_output = "rownames",
                                        input_species = "	
Crab-eating macaque",
                                        output_species = "Macaque",
                                        non121_strategy = "drop_both_species",
                                        method = "gprofiler")

## Get metadata from original object + create new object with orthologs only
write.csv(cyno_merge@meta.data, "metadata.csv")
meta <- read.csv("metadata.csv")
new_obj <- CreateAssayObject(gene_df)
cyno_Ortholog <- CreateSeuratObject(new_obj, meta.data = meta)

---------------------------------------------
## Normalize data
cyno_Ortholog <- NormalizeData(cyno_Ortholog, verbose = FALSE)
cyno_Ortholog <- FindVariableFeatures(cyno_Ortholog, selection.method = "vst", nfeatures = 2000)
cyno_Ortholog <- ScaleData(cyno_Ortholog, verbose = FALSE)
cyno_Ortholog <- RunPCA(cyno_Ortholog, npcs = 20, verbose = FALSE)
## ScaleData would usually also add in  vars.to.regress = c("S.Score","G2M.Score") - removed here

DimPlot(cyno_Ortholog, reduction = "umap", raster = FALSE, label = FALSE, group.by = c("dataset"))

---------------------------------------------
## Integrate data - Run Harmony
library(harmony)

cyno_Harmony <- RunHarmony(cyno_Ortholog, group.by.vars = "orig.ident", plot_convergence = TRUE)
cyno_Ortholog <- RunUMAP(cyno_Ortholog, reduction = "harmony", dims = 1:10)
cyno_Ortholog <- FindNeighbors(cyno_Ortholog, reduction = "harmony", dims = 1:10)
cyno_Ortholog <- FindClusters(object = cyno_Ortholog, resolution = 0.4, verbose = FALSE, random.seed = 1234)

saveRDS(cyno_Ortholog, "Cyno_RhesusOrthologs.rds")

---------------------------------------------
# Check feature plots
FeaturePlot(cyno_Ortholog, 
            raster = FALSE,
            cols = c("#D29343","#1d1147"),
            features = c("DDX4","DAZL","NANOS3","SYCE3","SOX17","TFAP2C","PRDM1","GATM","AMHR2","KRT19","NR2F2","TCF21")) &NoAxes()

ggsave("cyno_OrthologTestPanel.pdf", units="in", width = 12, height = 8, device='pdf', dpi=300)
ggsave("cyno_OrthologTestPanel.tiff", units="in", width = 12, height = 8, device='tiff', dpi=300)

# Call clusters - dotplot
# Define genes - used to call clusters in rhesus dotplot (Fig S4B)
whole_genes = c("SOX17","TFAP2C","DDX4","SYCP2","ZP3","KDR","PECAM1","VWF","UPK3B","KRT19","AMHR2","GATM","CD14","CD68","ACTA2","CNN1","NR2F1","NR2F2","TCF21","PDGFRA","SULT1E1")

# Plot
DefaultAssay(cyno_Ortholog) <- "RNA"
DotPlot(object = cyno_Ortholog, 
        features = whole_genes,
        dot.min = 0.01, 
        col.min = 0.0,
        cols = c("#D29343","#1d1147")) + RotatedAxis() + coord_flip() + xlab(' ')+ ylab(' ') +
  theme(legend.position = "right") +
  theme(legend.key.size = unit(0.5,"cm")) +
  scale_x_discrete(limits=rev)

ggsave("cyno_OrthologDotClusters.pdf", units="in", width = 6, height = 9, device='pdf', dpi=300)
ggsave("cyno_OrthologDotClusters.tiff", units="in", width = 6, height = 9, device='tiff', dpi=300)

# 13 clusters identified from FindClusters - celltype called as follows according to gene expression:
# Germline primordial: cluster 12, 5
# Germline meiotic: 2
# Granulosa: 0
# Stromal: 1, 3, 7
# Endothelial: 9
# Epithelial: 6
# Immune: 11
# Undefined: 4, 8, 10

# Organise clusters
Idents(cyno_Ortholog) <- "seurat_clusters"
table(Idents(cyno_Ortholog))
my_levels <- c(12,5,2,0,1,3,7,9,6,11,4,8,10)
Idents(cyno_Ortholog) <- factor(Idents(cyno_Ortholog), levels= my_levels3)
cyno_Ortholog[["cluster"]] <- Idents(cyno_Ortholog)

# Add celltype column
Idents(cyno_Ortholog) <- "seurat_clusters"
cyno_Ortholog <- RenameIdents(cyno_Ortholog, `12` = "Germline: primordial")
# Repeat for each cluster as detailed above, then save as column 
cyno_Ortholog[["celltype"]] <- Idents(cyno_Ortholog)

# Add dataset column with new nomenclature (e.g. W4_1(Sasaki))
Idents(cyno_Ortholog) <- "orig.ident"
table(Idents(cyno_Ortholog))
cyno_Ortholog <- RenameIdents(cyno_Ortholog, `Sasaki_E24` = "W4_1(Sasaki)")
# Repeat for each dataset as detailed above, then save as column 
cyno_Ortholog[["dataset"]] <- Idents(cyno_Ortholog)

saveRDS(cyno_Ortholog, "Cyno_RhesusOrthologs.rds")

---------------------------------------------
# Split datasets into embryonic (W4 - W6) and fetal (W8 - W18) granulosa
# Subset to select clusterX cells - W4 - W6
Idents(cyno_Ortholog) <- "dataset"
table(Idents(cyno_Ortholog))
sub_T1 <- subset(x = cyno_Ortholog, idents = c("W4_1(Sasaki)", "W4_2(Sasaki)","W5_1(Sasaki)","W5_2(Chen)","W5_3(Sasaki)","W6(Chen)"))

saveRDS(sub_T1, "Cyno_T1.rds")

# Subset to select clusterX cells - W8 - W18 Granulosa
sub_T2 <- subset(x = cyno_Ortholog, idents = c("W8_1(Mizuta)","W8_2(Mizuta)","W8_3(Chen)","W10(Mizuta)","W12(Mizuta)","W12(Zhao)","W16_1(Mizuta)","W16_2(Mizuta)","W17(Zhao)","W18(Mizuta)"))

Idents(sub_T2) <- "celltype"
sub_T2gran <- subset(x= sub_T2, idents = c("Granulosa"))
table(Idents(sub_T2gran))
#New table to confirm idents removed 

saveRDS(sub_T2gran, "Cyno_T2Gran.rds")

# Once subset, rerun standard clustering workflow to define new clusters
sub_T2gran <- FindVariableFeatures(sub_T2gran, selection.method = "vst", nfeatures = 2000)
sub_T2gran <- ScaleData(sub_T2gran, verbose = FALSE)
sub_T2gran <- RunPCA(sub_T2gran, npcs = 20, verbose = FALSE)
sub_T2gran <- RunUMAP(sub_T2gran, reduction = "harmony", dims = 1:10)
sub_T2gran <- FindNeighbors(sub_T2gran, reduction = "harmony", dims = 1:10)
sub_T2gran <- FindClusters(object = sub_T2gran, resolution = 0.2, verbose = FALSE, random.seed = 1234)

saveRDS(sub_T2gran, "Cyno_T2Gran.rds")
# Cluster resolution: sub_T1 = 0.2 , sub_T2Gran = 0.2

# For granulosa subset, cluster identity was defined based on gene expression: cluster 2 = Epithelial (KRT19, UPK3B),  cluster 1 = PG1 (NR2F2, PCP4, KDR), cluster 0 = PG2 (ITGA6, LHX2, LHX9, TOX3), cluster 3 = similar to rhesus cluster g8 (RDH10, HEYL, ETS1, TAGLN)
# Arrange clusters
my_levels2 <- c(2,0,1,3)
Idents(sub_T2gran) <- factor(Idents(sub_T2gran), levels= my_levels3)
sub_T2gran[["cluster"]] <- Idents(sub_T2gran)

saveRDS(sub_T2gran, "Cyno_T2Gran.rds")

-----------------------------
# Visualise data
# Cyno_T1 dataset UMAP with colour scale = viridis_rocket (from viridis)
DimPlot(sub_T1, pt.size = 0.2, group.by = "dataset", raster = FALSE, label = FALSE) +
  scale_color_viridis_d(option = "rocket", direction = -1)

ggsave("cyno_T1OrthologsDatasets.tiff", units="in", width = 12, height = 8, device='tiff', dpi=300)
ggsave("cyno_T1OrthologsDatasets.pdf", units="in", width = 12, height = 8, device='pdf', dpi=300

# Cyno_T1 celltype UMAP with manual colour values from batlow (khroma) palette
DimPlot(sub_T1, pt.size = 0.2, group.by = "celltype", raster = FALSE, label = TRUE) +
  scale_color_manual(values = c("#9c8831","#165a7f","#1e2256","#d19240","#edcce1","#f7b5bb","#aabc59","#F8A17B"))

# Cyno_T2 granulosa UMAP with batlow
DimPlot(cyno_Ortholog, pt.size = 0.2, group.by = "celltype", raster = FALSE, label = TRUE) +
  scale_color_batlow(discrete = TRUE, reverse = TRUE)

# Dotplots for Cyno_T1 and Cyno_T2 granulosa - adjust gene list as necessary
# Define genes
t1_genes <- c("DDX4", "ITGA6", "NR2F2", "LHX9", "NR5A1", "NR0B1", "EMX2", "KITLG", "AMHR2", "GATM", "RGS5", "DUOX2", "MFGE8", "MBNL3", "IFI6", "SERPINE2", "CYP17A1", "INHBA","GATA2", "PAX8","OSR2","OGN","LUM")
# Plot
DefaultAssay(sub_T1) <- "RNA"
DotPlot(object = sub_T1, 
        features = t1_genes,
        dot.min = 0.01, 
        col.min = 0.0,
        cols = c("#D29343","#1d1147")) + RotatedAxis() + coord_flip() + xlab(' ')+ ylab(' ') +
  theme(legend.position = "right") +
  theme(legend.key.size = unit(0.5,"cm")) +
  scale_x_discrete(limits=rev)

ggsave("cyno_T1OrthologDotClusters.pdf", units="in", width = 6, height = 9, device='pdf', dpi=300)
ggsave("cyno_T1OrthologDotClusters.tiff", units="in", width = 6, height = 9, device='tiff', dpi=300)
     