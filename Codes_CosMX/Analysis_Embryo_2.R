library(Seurat)
library(tibble)
library(progressr)
set.seed(1234)

nano.obj <- LoadNanostring(data.dir = "~/Desktop/CoSMX/Embryo_2", fov = "Embryo_2")

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
nano.obj <- RunPCA(nano.obj, dims = 1:30, verbose = FALSE)
nano.obj <- RunUMAP(nano.obj, dims = 1:30, verbose = FALSE)
nano.obj <- FindNeighbors(nano.obj, dims = 1:30, verbose = FALSE)
nano.obj <- FindClusters(nano.obj, verbose = FALSE, resolution = 5)
table(nano.obj@active.ident)
DimPlot(nano.obj, raster=FALSE)
ImageDimPlot(nano.obj, fov = "Embryo_2", axes = TRUE, cols = "glasbey")
AllMarkers <- FindAllMarkers(nano.obj, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/Embry_2_DEGs.csv")
write.csv(nano.obj@active.ident, "~/Desktop/Embryo_2_cell_type.txt")
AllMarkers <- FindAllMarkers(nano.obj, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/Embry_2_After_Annotation_DEGs.csv")

#Subset_1
Gonad_1 <- subset(nano.obj, (fov=="5")| (fov=="6")| (fov=="1")| (fov=="2")| (fov=="3")| (fov=="4")| (fov=="7")| (fov=="8")| (fov=="9"))
Gonad_1 <- SCTransform(Gonad_1, assay = "Nanostring", vst.flavor = "v2", clip.range = c(-10, 10), verbose = FALSE)
Gonad_1 <- RunPCA(Gonad_1, verbose = FALSE, npcs=100)
stdev <- Gonad_1@reductions$pca@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)

Gonad_1 <- RunUMAP(Gonad_1, dims = 1:PCNum, verbose = FALSE)
Gonad_1 <- FindNeighbors(Gonad_1, dims = 1:PCNum, verbose = FALSE)
Gonad_1 <- FindClusters(Gonad_1, verbose = FALSE, resolution = 3)
DimPlot(Gonad_1, raster=FALSE)
Idents(object = Gonad_1) <- "predicted.celltype"
ImageDimPlot(Gonad_1, fov = "Embryo_2", axes = TRUE, cols = "glasbey")
ggsave("~/Desktop/Gonad_1_Image_DimPlot.tiff", device = grDevices::tiff)

AllMarkers <- FindAllMarkers(Gonad_1, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/Gonad_1_DEGs.csv")

write.csv(Gonad_1@active.ident, "~/Desktop/Gonad_1_cell_type.txt")
Gonad_1@meta.data[["predicted.celltype"]] <- Gonad_1_cell_type$V2

batlowW <- c("#0C325E", "#1C5A62", "#FDB3B2", "#F8A483", "#F9CCF9", "#2E6A57", "#CE9243", "#77853A")

ImageDimPlot(Gonad_1, axes = TRUE, fov = "Embryo_2", group.by = "predicted.celltype") + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Gonad_1_annotated.tiff", device = grDevices::tiff)

ImageFeaturePlot(Gonad_1, fov = "Embryo_2",  features = c("NANOG", "KITLG", "KRT18", "KRT19", "NR2F2", "PDGFRA", "POU5F1", "PTGDS"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/CosMX_results_EMBRYO_2/Gonad_1_genes.tiff", device = grDevices::tiff)

ImageFeaturePlot(Gonad_1, fov = "Embryo_2",  features = c("MYL4", "MYL7", "OAS3", "RNF34", "CAV1", "MKI67"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/Gonad_1_genes_2.tiff", device = grDevices::tiff)

Idents(object = Gonad_1) <- "predicted.celltype"
AllMarkers <- FindAllMarkers(Gonad_1, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/After_annotation_Gonad_1_DEGs.csv")

###ZOOM OVARY#####
basal.crop <- Crop(Gonad_1[["Embryo_2"]], coords = "tissue", x = c(11000, 19000), y = c(110000, 123200))
Gonad_1[["zoom1"]] <- basal.crop
DefaultBoundary(Gonad_1[["zoom1"]]) <- "segmentation"

batlowW <- c("#77853A", "#F9CCF9", "#1C5A62", "#F8A483")

ImageDimPlot(Gonad_1, fov = "zoom1", coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Ovary_Zoom_Ovary_1.tiff", device = grDevices::tiff)


batlowW <- c("#1C5A62", "#77853A", "#F8A483", "#F9CCF9")
ImageDimPlot(Gonad_1, fov = "zoom1", group.by = "predicted.celltype", alpha = 0.3, molecules = c("RGS5", "IFI6", "KRT18", "KRT19"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Genes_Granulosa_Ovary_1.tiff", device = grDevices::tiff)

ImageDimPlot(Gonad_1, fov = "zoom1", group.by = "predicted.celltype", alpha = 0.3, molecules = c("NANOG", "POU5F1", "ITGA6"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Genes_Germline_Ovary_1.tiff", device = grDevices::tiff)

ImageDimPlot(Gonad_1, fov = "zoom1", group.by = "predicted.celltype", alpha = 0.3, molecules = c("NR2F2", "PDGFRA", "CYP1B1", "EFNB1"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Genes_Stroma_Ovary_1.tiff", device = grDevices::tiff)

################


#Subset_2
Gonad_2 <- subset(nano.obj, (fov=="10")| (fov=="11")| (fov=="12")| (fov=="13")| (fov=="14")| (fov=="15")| (fov=="16")| (fov=="17"))
Gonad_2 <- SCTransform(Gonad_2, assay = "Nanostring", vst.flavor = "v2", clip.range = c(-10, 10), verbose = FALSE)
Gonad_2 <- RunPCA(Gonad_2, verbose = FALSE, npcs = 100)

stdev <- Gonad_2@reductions$pca@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)

Gonad_2 <- RunUMAP(Gonad_2, dims = 1:PCNum, verbose = FALSE)
Gonad_2 <- FindNeighbors(Gonad_2, dims = 1:PCNum, verbose = FALSE)
Gonad_2 <- FindClusters(Gonad_2, verbose = FALSE, resolution = 2)
DimPlot(Gonad_2, raster=FALSE)
Idents(object = Gonad_2) <- "predicted.celltype"
ImageDimPlot(Gonad_2, fov = "Embryo_2", axes = TRUE, cols = "glasbey")
ggsave("~/Desktop/Gonad_2_Image_DimPlot.tiff", device = grDevices::tiff)

AllMarkers <- FindAllMarkers(Gonad_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/Gonad_2_DEGs.csv")

write.csv(Gonad_2@active.ident, "~/Desktop/Gonad_2_cell_type.txt")
Gonad_2@meta.data[["predicted.celltype"]] <- Gonad_2_cell_type$V2

batlowW <- c("#77853A",  "#1C5A62", "#F8A483", "#F9CCF9", "#F39E71", "#CE9243", "#A08A2D", "#CE9243", "#FDB3B2", "#E39855", "#557447", "#456E51", "#356859", "#275F5F", "#1D5460", "#144A60","#113E5F", "#0C325E")

batlowW <- c("#77853A", "#FDB3B2", "#1C5A62", "#0C325E", "#F8A483", "#CE9243", "#F9CCF9", "#2E6A57")

ImageDimPlot(Gonad_2, axes = TRUE, fov = "Embryo_2", group.by = "predicted.celltype") + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Gonad_2_annotated.tiff", device = grDevices::tiff)

ImageFeaturePlot(Gonad_2, fov = "Embryo_2",  features = c("NANOG", "KITLG", "KRT18", "KRT19", "NR2F2", "PDGFRA", "POU5F1", "PTGDS"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/CosMX_results_EMBRYO_2/Gonad_2_genes.tiff", device = grDevices::tiff)

ImageFeaturePlot(Gonad_2, fov = "Embryo_2",  features = c("MYL4", "MYL7", "OAS3", "RNF34", "CAV1", "MKI67"), max.cutoff = "q95") + NoLegend()
ggsave("~/Desktop/Gonad_2_genes_2.tiff", device = grDevices::tiff)

Idents(object = Gonad_2) <- "predicted.celltype"
AllMarkers <- FindAllMarkers(Gonad_2, only.pos = TRUE, logfc.threshold = 0)
write.csv(AllMarkers, "~/Desktop/CosMX_results_EMBRYO_2/After_annotation_Gonad_2_DEGs.csv")

###ZOOM OVARY#####
basal.crop <- Crop(Gonad_2[["Embryo_2"]], coords = "tissue", x = c(12500, 21000), y = c(7000, 16000))
Gonad_2[["zoom1"]] <- basal.crop
DefaultBoundary(Gonad_2[["zoom1"]]) <- "segmentation"
batlowW <- c("#77853A", "#F39E71", "#F9CCF9", "#1C5A62")

ImageDimPlot(Gonad_2, fov = "zoom1",  coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Ovary_Zoom_Ovary_2.tiff", device = grDevices::tiff)

batlowW <- c("#77853A", "#1C5A62", "#F39E71", "#F9CCF9")
ImageDimPlot(Gonad_2, fov = "zoom1", group.by = "predicted.celltype", alpha = 0.3, molecules = c("RGS5", "IFI6", "KRT18", "KRT19"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Genes_Granulosa_Ovary_2.tiff", device = grDevices::tiff)

ImageDimPlot(Gonad_2, fov = "zoom1", group.by = "predicted.celltype", alpha = 0.3, molecules = c("NANOG", "POU5F1", "ITGA6"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Genes_Germline_Ovary_2.tiff", device = grDevices::tiff)

ImageDimPlot(Gonad_2, fov = "zoom1", group.by = "predicted.celltype", alpha = 0.3, molecules = c("NR2F2", "PDGFRA", "CYP1B1", "EFNB1"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE) + scale_fill_manual(values = batlowW)
ggsave("~/Desktop/Genes_Stroma_Ovary_2.tiff", device = grDevices::tiff)



################
