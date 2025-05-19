library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

localdir <- "/Users/Macaca_HD/emb017/"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# Setting default assay to 8um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"

#Remove bins out of tissue
object@meta.data[["space_clusters"]] <- Only_tissue$Good
Idents(object) <- "space_clusters"
object <- subset(x = object, idents = c(0), invert = TRUE)

vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0, raster=FALSE) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")

vln.plot | count.plot

SpatialDimPlot(object, label = T, repel = T, label.size = 4)


# normalize 8um bins
object <- NormalizeData(object)

object <- FindVariableFeatures(object)
object <- ScaleData(object)
# we select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F, raster=FALSE) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2

SpatialDimPlot(object, label = T, repel = T, label.size = 4)


library(spacexr)
DefaultAssay(object) <- "sketch"

#import scRNA-seq
Idents(OvaryONPRC017) <- "Newcelltype"

#Remove cell type with less than 5 cells
OvaryONPRC017 <- subset(x = OvaryONPRC017, idents = c("Germline: primordial"), invert = TRUE)

counts <- OvaryONPRC017[["RNA"]]$counts
cluster <- as.factor(OvaryONPRC017$Newcelltype)
nUMI <- OvaryONPRC017$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD objecterence object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- object[["sketch"]]$counts
macaca_cells_hd <- colnames(object[["sketch"]])
coords <- GetTissueCoordinates(object)[macaca_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28, CELL_MIN_INSTANCE = 10, counts_MIN= 3, UMI_min = 6, gene_cutoff = 0, gene_cutoff_reg = 0, CONFIDENCE_THRESHOLD = 3, fc_cutoff_reg = 0.5)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
save(RCTD, file="/Users/RCTD.RData")
save(object, file="/Users/Emb_17_sc_predicted.RData")
# add results back to Seurat object
object <- AddMetaData(object, metadata = RCTD@results$results_df)

object$first_type <- as.character(object$first_type)
# project RCTD labels from sketched cortical cells to all cortical cells

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)


DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "full_first_type"

# now we can spatially map the location of any scRNA-seq cell type
cells <- CellsByIdentities(object)
cols <- c("Stromal" = "#77853A", " meiotic_Germline"= "#1C5A62", " granulosa_Somatic"= "#F9CCF9", "Endothelial"= "#F8A483", "Epithelial"= "#FDBACA", "Immune"="#FDB3B2", "Smooth muscle" ="#F29D6C")
p <- SpatialDimPlot(object, cols = cols)
p

####Markers####
markers <- FindAllMarkers(object, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_2 <- ScaleData(object, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_2, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p


#High_resolution cell types
library(spacexr)
DefaultAssay(object) <- "sketch"

#import scRNA-seq
OvaryONPRC017@meta.data[["High_resolution_cell_type"]] <- NEW_Ovary_cell__Meiotic$New_cells
Idents(OvaryONPRC017) <- "High_resolution_cell_type"

OvaryONPRC017 <- subset(x = OvaryONPRC017, idents = c("No"), invert = TRUE)
object_meiotic <- subset(x = object, idents = c(" meiotic_Germline"))


DefaultAssay(object_meiotic) <- "Spatial.008um"

counts <- OvaryONPRC017[["RNA"]]$counts
cluster <- as.factor(OvaryONPRC017$High_resolution_cell_type)
nUMI <- OvaryONPRC017$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD object_meioticerence object_meiotic
reference <- Reference(counts, cluster, nUMI)

counts_hd <- object_meiotic[["Spatial.008um"]]$counts
macaca_cells_hd <- colnames(object_meiotic[["Spatial.008um"]])
coords <- GetTissueCoordinates(object_meiotic)[macaca_cells_hd, 1:2]

# create the RCTD query object_meiotic
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28, CELL_MIN_INSTANCE = 10, counts_MIN= 3, UMI_min = 6, gene_cutoff = 0, gene_cutoff_reg = 0, CONFIDENCE_THRESHOLD = 3, fc_cutoff_reg = 0.5)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
save(RCTD, file="/Users/RCTD_High_resolution_cell_type.RData")
save(object_meiotic, file="/Users/Emb_17_sc_predicted_High_resolution_cell_type.RData")
# add results back to Seurat object_meiotic
object_meiotic <- AddMetaData(object_meiotic, metadata = RCTD@results$results_df)

object_meiotic$first_type <- as.character(object_meiotic$first_type)

DefaultAssay(object_meiotic) <- "Spatial.008um"
Idents(object_meiotic) <- "first_type"

p <- SpatialDimPlot(object_meiotic)
p

#import scRNA-seq
OvaryONPRC017@meta.data[["High_resolution_cell_type"]] <- NEW_Ovary_cell_granulosa$New_cells
Idents(OvaryONPRC017) <- "High_resolution_cell_type"

OvaryONPRC017 <- subset(x = OvaryONPRC017, idents = c("No", "Theca-like t1"), invert = TRUE)
object_granulosa <- subset(x = object, idents = c(" granulosa_Somatic"))

DefaultAssay(object_granulosa) <- "Spatial.008um"

counts <- OvaryONPRC017[["RNA"]]$counts
cluster <- as.factor(OvaryONPRC017$High_resolution_cell_type)
nUMI <- OvaryONPRC017$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD object_granulosaerence object_granulosa
reference <- Reference(counts, cluster, nUMI)

counts_hd <- object_granulosa[["Spatial.008um"]]$counts
macaca_cells_hd <- colnames(object_granulosa[["Spatial.008um"]])
coords <- GetTissueCoordinates(object_granulosa)[macaca_cells_hd, 1:2]

# create the RCTD query object_granulosa
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28, CELL_MIN_INSTANCE = 10, counts_MIN= 3, UMI_min = 6, gene_cutoff = 0, gene_cutoff_reg = 0, CONFIDENCE_THRESHOLD = 3, fc_cutoff_reg = 0.5)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
save(RCTD, file="/Users/RCTD_High_resolution_cell_type.RData")
save(object_granulosa, file="/Users/Emb_17_sc_predicted_High_resolution_cell_type.RData")
# add results back to Seurat object_granulosa
object_granulosa <- AddMetaData(object_granulosa, metadata = RCTD@results$results_df)

object_granulosa$first_type_granulosa <- as.character(object_granulosa$first_type)

DefaultAssay(object_granulosa) <- "Spatial.008um"
Idents(object_granulosa) <- "first_type_granulosa"

p <- SpatialDimPlot(object_granulosa)
p