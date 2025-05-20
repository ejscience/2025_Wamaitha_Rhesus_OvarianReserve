# Analysis of rhesus datasets generated for this study -------------------------------
library(Seurat)
library(tidyverse)
library(DropletUtils)
library(BiocManager)
library(patchwork)
library(viridis)
library(khroma)
library(dittoSeq)
library(wesanderson)
library(MAST)
library(pheatmap)


--------------------------
# QC and processing
## Load merged seurat object 
ovary_all <- readRDS("~/Box/Genomics/10x Oregon analysis/Ovary/all18ov_v3.seurat.rds")

## QC - cell cycle scoring, then subset on nFeature, nCount, percent.mito (percent.mt)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ovary_all <- CellCycleScoring(ovary_all, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## Check and save QC plots
VlnPlot(ovary_whole, group.by = "animalID", features = c("nFeature_RNA","nCount_RNA","p.mito","percent.ribo"), pt.size = 0.1, ncol = 2)
ggsave("ovaryWhole_QCplots.pdf", units="in", width = 8, height = 8, device='pdf', dpi=300)
## Subset
ovary_all <- subset(ovary_all, subset = nFeature_RNA > 500)
ovary_all <- subset(ovary_all, subset = nCount_RNA > 1000)
ovary_all <- subset(ovary_all, subset = percent.mt < 15)
## ONPRC017_A excluded due to low QC metrics
### Idents(ovary_whole) <- "DatasetName"
### table(Idents(ovary_whole))
### ovary_whole <- subset(x = ovary_whole, idents = ONPRCO17_XX_D130_A_S9, invert = TRUE)

## Split the dataset into a list of Seurat objects
ovary.list <- SplitObject(ovary_all, split.by = "orig.ident")

## Normalize and identify variable features for each dataset independently
ovary.list <- lapply(X = ovary.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
## Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ovary.list)

## Normalize and identify variable features for each dataset independently
ovary_all <- NormalizeData(ovary_all, verbose = FALSE)
ovary_all <- FindVariableFeatures(ovary_all, selection.method = "vst", nfeatures = 2000)

## Seurat v4 integration
ovary.anchors <- FindIntegrationAnchors(object.list = ovary.list, anchor.features = features,     reduction = "rpca")
ovary_integrated <- IntegrateData(anchorset = ovary.anchors)

## Integrated analysis
DefaultAssay(ovary_integrated) <- "integrated"
ovary_integrated <- withCallingHandlers(expr, warning = function(w) if (inherits(w,     classes)) tryInvokeRestart("muffleWarning"), reduction = "pca", dims = 1:30)
ovary_integrated <- ScaleData(ovary_integrated, vars.to.regress = c("p.mito","S.Score","G2M.Score"), verbose = FALSE)
ovary_integrated <- RunPCA(ovary_integrated, npcs = 30, verbose = FALSE)
ovary_integrated <- RunUMAP(ovary_integrated, reduction = "pca", dims = 1:30)
ovary_integrated <- FindNeighbors(ovary_integrated, reduction = "pca", dims = 1:30)
ovary_integrated <- FindClusters(object = ovary_integrated, resolution = 0.15, verbose = FALSE, random.seed = 1234)

saveRDS(ovary_integrated, "ovary_integrated.rds")

## Add dataset column to group datasets by timepoint 
Idents(ovary_integrated) <- "animalID"
ovary_integrated <- RenameIdents(ovary_integrated, `ONPRC001` = "W19")
## Repeat for each sample to create W8, W15 and W19 groups, then save as column 
ovary_integrated[["group"]] <- Idents(ovary_integrated)

## Specify sample group order
group_order <- c("W8","W15","W19")
ovary_integrated$group <- factor(ovary_integrated$group, levels = group_order)


saveRDS(ovary_integrated, "ovary_integrated.rds")

---------------------------------
# Cluster analysis
## 16 clusters identified from FindClusters - celltype called as follows according to gene expression:
## Germline primordial: cluster a12
## Germline meiotic: a6
## Granulosa: a1, a3, a4
## Smooth muscle: a11
## Stromal: a0, a2, a5. a10
## Endothelial: a7, a 15
## Epithelial: a14
## Immune: a13
## Undefined: a8, a9

# Visualise UMAP clusters - viridis or library(khroma) Crameri colour scheme: https://www.nature.com/articles/s41467-020-19160-7

DimPlot(ovary_integrated, pt.size = 0.1, group.by = "animalID", raster = FALSE, label = FALSE) +
  scale_color_manual(values = c("#f9cb35","#fec488","#fc8961","#e75263","#b73779","#832681","#51127c","#1d1147","#000004" ), 
                     aesthetics = "colour") &NoLegend()
ggsave("ovaryWhole_animal.pdf", units="in", width = 6, height = 9, device='pdf', dpi=300)
## Colour scale = viridis_magma with inferno first step

DimPlot(ovary_integrated, pt.size = 0.2, group.by = "seurat_clusters", raster = FALSE, label = TRUE) +
  scale_colour_manual(values = c("#FACCFA","#FDBCCB","#FDB7BC","#F8A17B","#FCA993","#FDB1AB","#ED9A62","#D29343","#9D892B","#BE9035","#FAA587","#687B3E","#2D675D","#103F60","#011959","#DE964F")) &NoLegend()

"#FACCFA"- stroma1_0
"#FDB7BC"- stroma2_2
"#FDB1AB"- stroma3_5
"#FAA587"- stroma4_10
"#FDBCCB"- granulosa1_1
"#FCA993" - granulosa3_4
"#9D892B" - undefined1_8
"#BE9035" - undefined2_9
"#D29343" - endothelial_7
"#DE964F" - endothelial2_15
"#F8A17B" - granulosa2_3
"#ED9A62" - germline:meiotic_6
"#2D675D" - germline:primordial/all_12
"#687B3E" - smoooth:muscle_11
"#103F60" - immune_13
"#011959" - epithelial_14

ggsave("ovaryWhole_clusters.pdf", units="in", width = 6, height = 9, device='pdf', dpi=300)
## To split by timepoint e.g. Figure S4C  add 'split.by = "group"' command in DimPlot syntax

## Add celltype column to group datasets by timepoint 
Idents(ovary_integrated) <- "seurat_clusters"
ovary_integrated <- RenameIdents(ovary_integrated, `12` = "Germline: primordial")
## Repeat for each cluster, then save as column 
ovary_integrated[["celltype"]] <- Idents(ovary_integrated)

saveRDS(ovary_integrated, "ovary_integrated.rds")

# Dotplot cluster gene expression
## Define genes
whole_genes = c("SOX17","TFAP2C","DDX4","SYCP2","ZP3","KDR","PECAM1","VWF","UPK3B","KRT19","AMHR2","GATM","CD14","CD68","ACTA2","CNN1","NR2F1","NR2F2","TCF21","PDGFRA","SULT1E1")

## Plot
DefaultAssay(ovary_integrated) <- "RNA"
DotPlot(object = ovary_integrated, 
        features = whole_genes,
        dot.min = 0.01, 
        col.min = 0.0,
        cols = c("#D29343","#1d1147")) + RotatedAxis() + coord_flip() + xlab(' ')+ ylab(' ') +
  theme(legend.position = "right") +
  theme(legend.key.size = unit(0.5,"cm")) +
  scale_x_discrete(limits=rev)

ggsave("ovary_wholeDotClusters.pdf", units="in", width = 6, height = 9, device='pdf', dpi=300)

# Bar graph cluster proportions - dittoseq
dittoBarPlot(
  object = ovary_integrated,
  var = "celltype",
  group.by = "group",
  x.reorder = c(7,8,9,2,3,4,1,5,6),
  color.panel = c("#D29343","#011959","#F8A17B","#3C6D56","#FDB7BC","#103F60","#687B3E","#FACCFA","#9D892B"))
ggsave("ovaryWhole_clusterBarAnimal.pdf", units="in", width = 4, height = 4, device='pdf', dpi=300)


# Violin plots
VlnPlot(ovary_integrated, features = c("ESR1","ESR2","GPER1","AR"), stack = TRUE, flip = TRUE, raster = FALSE, split.by = "group") +
  theme(legend.position = "none") +
  ggtitle(" ")

ggsave("ovarywhole_hormoneVln.pdf", units="in", width = 8, height = 12, device='pdf', dpi=300)


---------------------------------------------
# Marker analysis
## Find Conserved Markers
get_conserved <- function(cluster){
  FindConservedMarkers(ovary_integrated,
                       ident.1 = cluster,
                       grouping.var = "animalID",
                       only.pos = TRUE,
                       #thresholds used by Lina
                       logfc.threshold = 0,min.pct = 0) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

##set params
DefaultAssay(ovary_integrated) <- "RNA"
Idents(ovary_integrated) <-  "integrated_snn_res.0.15"

conserved_markers <- map_dfr(0:15, get_conserved) #map across clusters 0-15
saveRDS(conserved_markers, "conserved_markers_res0.15.rds")

write_csv(conserved_markers %>% 
            arrange(cluster_id, desc(max_pval)),
          file = "markers_for_each_cluster_res.0.15.csv")

markers <- readRDS("conserved_markers_res0.15.rds")

#print scrollable list
markers %>% group_by(cluster_id) %>% 
  select(gene, max_pval, everything()) %>% slice_head(n = 20) 


##Find markers for every cluster compared to all remaining cells, report only the positives
library(MAST)

ovary.markers <- FindAllMarkers(object = ovary_integrated, only.pos = TRUE, test.use = "MAST")
ovaryMast.csv <- print(ovary.markers %>% group_by(cluster) %>% top_n(100, avg_log2fc))
write.csv (ovaryMast.csv, "ovaryWhole_mast.csv")


-----------------------------------
# Subset Analysis
## Subsample by group - W8
obj.list <- SplitObject(ovary_integrated, split.by = "group")
saveRDS(obj.list$W8, "W8_subset.rds")
##Once subset, rerun standard workflow (without integration step)
##W8 cluster resolution 0.15

## Subsample by cluster - e.g. germline cells
table(Idents(ovary_integrated))
sub_germ <- subset(x = ovary_integrated, idents = c("Germline: primordial", "Germline: meiotic"))
table(Idents(sub_germ))
##New table to confirm idents removed 
saveRDS(sub_germ, "Ovary_germline.rds")
##Identical process for granulosa subset, just switch idents in command
##Once subset, rerun standard workflow (without integration step)
##Germline cluster resolution 0.2 ; Granulosa resolution 0.3 ; Stroma resolution 0.15 ; Theca resolution 0.15 (subset from reclustered Granulosa subset)
  
## Check UMAP feature plots
FeaturePlot(ovary_W8, 
              raster = FALSE,
              cols = c("#D29343","#1d1147"),
              features = c("GATM","AMHR2","LHX9","RGS5","MFGE8","NR0B1")) &NoAxes()
ggsave("ovaryW8_CytAssistUp.pdf", units="in", width = 6, height = 9, device='pdf', dpi=300)


## Visualise UMAP clusters for theca subset Fig S7C - library (wesanderson) Royal2 colour palette 
DimPlot(sub_granTheca, pt.size = 2, raster = FALSE, label = TRUE) +
  scale_color_manual(values = wes_palette(3, name="Royal2")) & NoLegend()


---------------------------
# Pseudobulk analysis
## Make heatmap with pheatmap - tutorial at bioinformatics.ccr.cancer.gov/docs/data-visualization-with-r/Lesson5_intro_to_ggplot/#import-data
library(pheatmap)
mat <- read.csv("~/Box/Data/Genomics/10x Oregon analysis/Ovary_granulosa/pseudoBulk_granulosa/W8_pseudobulkDE_FDRadjVartop500.csv", header = TRUE, row.names = 1, check.names = FALSE)

p2 <- as.ggplot(pheatmap(mat, 
                         scale = "row", 
                         color = colorRampPalette(c("#D29343","#1d1147"))(25),
                         show_rownames = FALSE,
                         cellwidth = 12,
                         treeheight_col = 25,
                         treeheight_row = 25,
                         cutree_cols = 3, 
                         cutree_rows = 6))


saveRDS(pheatGran_top10pseudo, file="~/Box/Data/Genomics/10x Oregon analysis/Ovary_W8/pseudoBulk_W8/W8_pseudobulkDE_FDRadjVartop500.rds")

### Pull genes from k-means clusters:
hm <- (pheatmap(mat, 
                scale = "row", 
                color = colorRampPalette(c("#D29343","#1d1147"))(25),
                show_rownames = FALSE,
                cellwidth = 12,
                treeheight_col = 25,
                treeheight_row = 25,
                cutree_cols = 3, 
                cutree_rows = 6))

row_cluster = data.frame(cluster = cutree(hm$tree_row, k = 6))
write.csv(row_cluster, 'pheatmap_clusters.csv')
### add annotation_row = row_cluster to p2 code to get cluster labels

ggsave("pheatW8_pseudoCluster.tiff", units="in", width = 6, height = 9, device='tiff', dpi=300)