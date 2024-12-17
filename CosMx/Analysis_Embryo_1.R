library(Seurat)
library(tibble)
library(progressr)
library(tidyverse)
library(rgeos)
set.seed(1234)

nano.obj <- LoadNanostring(data.dir = "~/Desktop/CoSMX/Embryo_1", fov = "Embryo_1")

nano.obj@meta.data <- rownames_to_column(nano.obj@meta.data, var = 'rowname')
nano.obj@meta.data$fov <- nano.obj@meta.data$rowname
nano.obj@meta.data <- column_to_rownames(nano.obj@meta.data, var = 'rowname')
nano.obj@meta.data$fov <- gsub(".*_","", nano.obj@meta.data$fov)

# Visualize QC metrics as a violin plot
dim(nano.obj)
nano.obj[["percent.SysPrb"]] <- PercentageFeatureSet(nano.obj, pattern = "^SystemControl\\d{1,3}$")
nano.obj[["percent.NegPrb"]] <- PercentageFeatureSet(nano.obj, pattern = "^Negative\\d{1,3}$")
VlnPlot(nano.obj, features = c("nFeature_Nanostring", "nCount_Nanostring", "percent.NegPrb", "percent.SysPrb"), ncol = 4)
nano.obj <- subset(nano.obj, subset = nFeature_Nanostring > 20 & percent.NegPrb < 25 & percent.SysPrb < 25)
VlnPlot(nano.obj, features = c("nFeature_Nanostring", "nCount_Nanostring", "percent.NegPrb", "percent.SysPrb"), ncol = 4)
dim(nano.obj)
nano.obj <- nano.obj[!grepl("^SystemControl(\\d{1,3})$", rownames(nano.obj)), ]
nano.obj <- nano.obj[!grepl("^Negative(\\d{1,3})$", rownames(nano.obj)), ]
dim(nano.obj)


#Analysis
nano.obj <- SCTransform(nano.obj, assay = "Nanostring", vst.flavor = "v2", clip.range = c(-10, 10), verbose = FALSE)
nano.obj <- RunPCA(nano.obj, verbose = FALSE)
nano.obj <- RunUMAP(nano.obj, dims = 1:30, verbose = FALSE)
nano.obj <- FindNeighbors(nano.obj, dims = 1:30, verbose = FALSE)
nano.obj <- FindClusters(nano.obj, verbose = FALSE)
DimPlot(nano.obj, raster=FALSE)
ImageDimPlot(nano.obj, fov = "Embryo_1", axes = TRUE, cols = "glasbey")
AllMarkers <- FindAllMarkers(nano.obj, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/Nano_1_DEGs.csv")
write.csv(nano.obj@active.ident, "~/Desktop/nano.obj_cell_type.txt")
nano.obj@meta.data[["predicted.celltype"]] <- nano.obj_cell_type$V2
Idents(nano.obj) <- "predicted.celltype"


#Subset_1
Gonad_1 <- subset(nano.obj, (fov=="1")| (fov=="4")| (fov=="5")| (fov=="7")| (fov=="8")| (fov=="9")| (fov=="10"))
Gonad_1 <- SCTransform(Gonad_1, assay = "Nanostring", vst.flavor = "v2", clip.range = c(-10, 10), verbose = FALSE)
Gonad_1 <- RunPCA(Gonad_1, verbose = FALSE)
Gonad_1 <- RunUMAP(Gonad_1, dims = 1:30, verbose = FALSE)
Gonad_1 <- FindNeighbors(Gonad_1, dims = 1:30, verbose = FALSE)
Gonad_1 <- FindClusters(Gonad_1, verbose = FALSE)
DimPlot(Gonad_1, raster=FALSE)
ImageDimPlot(Gonad_1, axes = TRUE, fov = "Embryo_1", cols = "glasbey")
ggsave("~/Desktop/CosMX_results_EMBRYO_1/test.tiff", device = grDevices::tiff)
dev.off()
AllMarkers <- FindAllMarkers(Gonad_1, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/CosMX_results_EMBRYO_1/Gonad_1_DEGs.csv")

write.csv(Gonad_1@active.ident, "~/Desktop/Gonad_1_cell_type.txt")
Gonad_1@meta.data[["predicted.celltype"]] <- Gonad_1_cell_type$V2
Idents(Gonad_1) <- "predicted.celltype"

ImageDimPlot(Gonad_1, axes = TRUE, fov = "Embryo_1", cols = batlowW, group.by = "predicted.celltype")
ggsave("~/Desktop/Gonad_1_annotated.tiff", device = grDevices::tiff)

batlow <- color("batlow")
plot(batlow(256))

batlowW <- c("#77853A", "#1C5A62", "#F8A483", "#F9CCF9", "#F39E71", "#CE9243", "#A08A2D", "#CE9243", "#FDB3B2", "#E39855", "#557447", "#456E51", "#356859", "#275F5F", "#1D5460", "#144A60","#113E5F", "#0C325E")

ImageDimPlot(Gonad_1, axes = TRUE, fov = "Embryo_1", group.by = "predicted.celltype") + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Gonad_1_annotated.tiff", device = grDevices::tiff)

###ZOOM OVARY#####
basal.crop <- Crop(Gonad_1[["Embryo_1"]], coords = "tissue", x = c(900, 4600), y = c(107200, 113500))
Gonad_1[["zoom1"]] <- basal.crop
DefaultBoundary(Gonad_1[["zoom1"]]) <- "segmentation"

batlowW <- c("#77853A", "#F8A483","#1C5A62", "#F9CCF9", "#F39E71", "#CE9243", "#A08A2D", "#CE9243", "#FDB3B2", "#E39855", "#557447", "#456E51", "#356859", "#275F5F", "#1D5460", "#144A60","#113E5F", "#0C325E")
ImageDimPlot(Gonad_1, fov = "zoom1", cols = batlowW, coord.fixed = FALSE)
ggsave("~/Desktop/Ovary_Zoom.tiff", device = grDevices::tiff)

batlowW <- c("#1C5A62", "#F9CCF9","#77853A", "#F8A483", "#F39E71", "#CE9243", "#A08A2D", "#CE9243", "#FDB3B2", "#E39855", "#557447", "#456E51", "#356859", "#275F5F", "#1D5460", "#144A60","#113E5F", "#0C325E")
ImageDimPlot(Gonad_1, fov = "zoom1", cols = batlowW, group.by = "predicted.celltype", alpha = 0.5, molecules = c("RGS5", "IFI6", "KRT18", "KRT19"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)
ggsave("~/Desktop/Genes_Granulosa.tiff", device = grDevices::tiff)

ImageDimPlot(Gonad_1, fov = "zoom1", cols = batlowW, group.by = "predicted.celltype", alpha = 0.5, molecules = c("NANOG", "ITGA6", "POU5F1"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)
ggsave("~/Desktop/Genes_Germline.tiff", device = grDevices::tiff)

ImageDimPlot(Gonad_1, fov = "zoom1", cols = batlowW, group.by = "predicted.celltype", alpha = 0.5, molecules = c("NR2F2", "PDGFRA", "CYP1B1", "EFNB1"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)
ggsave("~/Desktop/Stroma_Germline.tiff", device = grDevices::tiff)
################

ImageFeaturePlot(Gonad_1, fov = "Embryo_1",  features = c("NANOG", "KITLG", "KRT18", "KRT19", "NR2F2", "PDGFRA", "POU5F1", "PTGDS"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/CosMX_results_EMBRYO_1/Gonad_1_genes.tiff", device = grDevices::tiff)

ImageFeaturePlot(Gonad_1, fov = "Embryo_1",  features = c("MYL4", "MYL7", "OAS3", "RNF34", "CAV1", "MKI67"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/Gonad_1_genes2.tiff", device = grDevices::tiff)

Idents(object = Gonad_1) <- "predicted.celltype"
AllMarkers <- FindAllMarkers(Gonad_1, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/CosMX_results_EMBRYO_1/After_annotation_Gonad_1_DEGs.csv")

#Subset_2
Gonad_2 <- subset(nano.obj, (fov=="11")| (fov=="12")| (fov=="13"))
Gonad_2 <- SCTransform(Gonad_2, assay = "Nanostring", vst.flavor = "v2", clip.range = c(-10, 10), verbose = FALSE)
Gonad_2 <- RunPCA(Gonad_2, verbose = FALSE)
Gonad_2 <- RunUMAP(Gonad_2, dims = 1:30, verbose = FALSE)
Gonad_2 <- FindNeighbors(Gonad_2, dims = 1:30, verbose = FALSE)
Gonad_2 <- FindClusters(Gonad_2, verbose = FALSE)
DimPlot(Gonad_2, raster=FALSE)
ImageDimPlot(Gonad_2, fov = "Embryo_1", axes = TRUE, cols = "glasbey")
ggsave("~/Desktop/CosMX_results_EMBRYO_1/Gonad_2_Image_DimPlot.tiff", device = grDevices::tiff)
AllMarkers <- FindAllMarkers(Gonad_2, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/CosMX_results_EMBRYO_1/Gonad_2_DEGs.csv")

write.csv(Gonad_2@active.ident, "~/Desktop/Gonad_2_cell_type.txt")
Gonad_2@meta.data[["predicted.celltype"]] <- Gonad_2_cell_type$V2

ImageDimPlot(Gonad_2, axes = TRUE, fov = "Embryo_1", cols = "glasbey", group.by = "predicted.celltype")
ggsave("~/Desktop/CosMX_results_EMBRYO_1/Gonad_2_annotated.tiff", device = grDevices::tiff)

ImageFeaturePlot(Gonad_2, fov = "Embryo_1",  features = c("NANOG", "KITLG", "KRT18", "KRT19", "NR2F2", "PDGFRA", "POU5F1", "PTGDS"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/CosMX_results_EMBRYO_1/Gonad_2_genes.tiff", device = grDevices::tiff)

Idents(object = Gonad_2) <- "predicted.celltype"
AllMarkers <- FindAllMarkers(Gonad_2, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/CosMX_results_EMBRYO_1/After_annotation_Gonad_2_DEGs.csv")

###ZOOM OVARY#####
basal.crop <- Crop(Gonad_2[["Embryo_1"]], coords = "tissue", x = c(11000, 15000), y = c(-1000, 3000))
Gonad_2[["zoom1"]] <- basal.crop
DefaultBoundary(Gonad_2[["zoom1"]]) <- "segmentation"

ImageDimPlot(Gonad_2, fov = "zoom1", cols = "polychrome", coord.fixed = FALSE)
ImageDimPlot(Gonad_2, fov = "zoom1", cols = "polychrome", group.by = "predicted.celltype", alpha = 0.3, molecules = c("RGS5", "IFI6"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)
############

#Subset_3
Gonad_3 <- subset(nano.obj, (fov=="2")| (fov=="6")| (fov=="3"))
Gonad_3 <- SCTransform(Gonad_3, assay = "Nanostring", vst.flavor = "v2", clip.range = c(-10, 10), verbose = FALSE)
Gonad_3 <- RunPCA(Gonad_3, verbose = FALSE)
Gonad_3 <- RunUMAP(Gonad_3, dims = 1:30, verbose = FALSE)
Gonad_3 <- FindNeighbors(Gonad_3, dims = 1:30, verbose = FALSE)
Gonad_3 <- FindClusters(Gonad_3, verbose = FALSE)
DimPlot(Gonad_3, raster=FALSE)
ImageDimPlot(Gonad_3, fov = "Embryo_1", axes = TRUE, cols = "glasbey")
ggsave("~/Desktop/CosMX_results_EMBRYO_1/Gonad_3_Image_DimPlot.tiff", device = grDevices::tiff)
AllMarkers <- FindAllMarkers(Gonad_3, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/CosMX_results_EMBRYO_1/Gonad_3_DEGs.csv")

write.csv(Gonad_3@active.ident, "~/Desktop/Gonad_3_cell_type.txt")
Gonad_3@meta.data[["predicted.celltype"]] <- Gonad_3_cell_type$V2

ImageDimPlot(Gonad_3, axes = TRUE, fov = "Embryo_1", cols = "glasbey", group.by = "predicted.celltype")
ggsave("~/Desktop/CosMX_results_EMBRYO_1/Gonad_3_annotated.tiff", device = grDevices::tiff)

ImageFeaturePlot(Gonad_3, fov = "Embryo_1",  features = c("NANOG", "KITLG", "KRT18", "KRT19", "NR2F2", "PDGFRA", "POU5F1", "PTGDS"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/CosMX_results_EMBRYO_1/Gonad_3_genes.tiff", device = grDevices::tiff)

Idents(object = Gonad_3) <- "predicted.celltype"
AllMarkers <- FindAllMarkers(Gonad_3, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/CosMX_results_EMBRYO_1/After_annotation_Gonad_3_DEGs.csv")

###ZOOM OVARY#####
basal.crop <- Crop(Gonad_3[["Embryo_1"]], coords = "tissue", x = c(14000, 18000), y = c(120500, 124000))
Gonad_3[["zoom1"]] <- basal.crop
DefaultBoundary(Gonad_3[["zoom1"]]) <- "segmentation"

ImageDimPlot(Gonad_3, fov = "zoom1", cols = "polychrome", coord.fixed = FALSE)
ImageDimPlot(Gonad_3, fov = "zoom1", cols = "polychrome", group.by = "predicted.celltype", alpha = 0.3, molecules = c("RGS5", "IFI6"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)
############

