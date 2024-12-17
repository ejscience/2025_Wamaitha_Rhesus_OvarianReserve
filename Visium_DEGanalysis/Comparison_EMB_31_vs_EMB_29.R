library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

#Load
emb_31 <- Load10X_Spatial('~/Desktop/New_Visium_Rm/emb031_outs/outs/',filename = "filtered_feature_bc_matrix.h5", slice = "slice1")

#QC
VlnPlot(emb_31, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(emb_31, features = c("nCount_Spatial", "nFeature_Spatial"))
emb_31 <- emb_31[, emb_31$nFeature_Spatial > 100 & emb_31$nCount_Spatial > 100]
SpatialFeaturePlot(emb_31, features = c("nCount_Spatial", "nFeature_Spatial"))
dim(emb_31)

#SCT trasnform
emb_31 <- SCTransform(emb_31, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)
emb_31 <- RunPCA(emb_31, assay = "SCT", verbose = FALSE)
emb_31 <- FindNeighbors(emb_31, reduction = "pca", dims = 1:30)
emb_31 <- FindClusters(emb_31, verbose = FALSE)
emb_31 <- RunUMAP(emb_31, reduction = "pca", dims = 1:30)

p1 <- DimPlot(emb_31, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(emb_31, label = TRUE, label.size = 1, pt.size.factor = 1)
p1 + p2

cells_7 <- c("TCATAGCTTCCAGTAC-1",
             "GCGGTGGCCGACGCAT-1",
             "CTCAATGACCATTAAG-1",
             "GGCGTGTAATTACTTA-1",
             "TAAGTGGCCGCTCACC-1",
             "ATCCGGATGACGTAAC-1",
             "GGTCCTGCGAGTAGGT-1",
             "AAGTTCCTATTGCTGG-1",
             "GAACATTCGACACGTA-1",
             "CGGTATAGGATGAGCG-1",
             "AACCACTGCCATAGCC-1",
             "CAAGCCGCCATAGAGG-1"
)



Idents(emb_31) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_31, cells_7) <- "Ovaries"


SpatialDimPlot(emb_31, pt.size.factor = 1)

####EMB_29####
#Load
emb_29 <- Load10X_Spatial('~/Desktop/NEW_Visium_Rm/emb029_outs/outs/',filename = "filtered_feature_bc_matrix.h5", slice = "slice3")
#QC
dim(emb_29)
VlnPlot(emb_29, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(emb_29, features = c("nCount_Spatial", "nFeature_Spatial"))
emb_29 <- emb_29[, emb_29$nFeature_Spatial > 100 & emb_29$nCount_Spatial > 100]
SpatialFeaturePlot(emb_29, features = c("nCount_Spatial", "nFeature_Spatial"))
dim(emb_29)

#SCT trasnform
emb_29 <- SCTransform(emb_29, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)
emb_29 <- RunPCA(emb_29, assay = "SCT", verbose = FALSE)
emb_29 <- FindNeighbors(emb_29, reduction = "pca", dims = 1:30)
emb_29 <- FindClusters(emb_29, verbose = FALSE)
emb_29 <- RunUMAP(emb_29, reduction = "pca", dims = 1:30)

p1 <- DimPlot(emb_29, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(emb_29, label = TRUE, label.size = 3, pt.size.factor = 1)
p1 + p2

cells_5 <- c(
  "AACCGCCAGACTACTT-1",
  "CGCAGATAGATGTTCT-1",
  "CTACCAGACCTCTTAG-1",
  "GACAATTGTATGCTTC-1",
  "TCGGCGCAGATAGATA-1",
  "TACTCCGAAGTAGAAT-1",
  "TAACAGGTCCATACCA-1",
  "TCCGAGTACGACTTAG-1",
  "GGTTCACCGACTCACG-1",
  "TAGTTCGAGTAGGAAG-1",
  "ATTAGTGCTTGGTCAT-1",
  "CGTAGTACGCATAATG-1",
  "TCATAGCTTCCAGTAC-1",
  "GGTCTTGAACCGGCCA-1",
  "TAAGTGGCCGCTCACC-1",
  "TCCGCCAGTGATTAAT-1",
  "TCAGGAGCAACGGTCG-1",
  "GCGGTGGCCGACGCAT-1",
  "GGCATTGCTAAGAATG-1",
  "TCTTACTCTGTATGCG-1",  
  "AATGGACATCCTACTC-1"
)

Idents(emb_29) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_29, cells_5) <- "Testis"

SpatialDimPlot(emb_29, pt.size.factor = 1)



emb_29@meta.data[["orig.ident"]] <- "emb_29"
emb_31@meta.data[["orig.ident"]] <- "emb_31"

#Integration
embryo_list <- list(Object_name1 = emb_31, 
                    Object_name2 = emb_29)

st.features1 = SelectIntegrationFeatures(embryo_list, nfeatures = 3000, verbose = FALSE)
embryo_list <- PrepSCTIntegration(object.list = embryo_list, anchor.features = st.features1)
combined.anchors1 <- FindIntegrationAnchors(object.list = embryo_list, normalization.method = "SCT",
                                            anchor.features = st.features1)
embryo.integrated <- IntegrateData(anchorset = combined.anchors1, normalization.method = "SCT")
embryo.integrated <- RunPCA(embryo.integrated, verbose = FALSE)

stdev <- embryo.integrated@reductions$pca@stdev
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

embryo.integrated <- FindNeighbors(embryo.integrated, dims = 1:PCNum)
embryo.integrated <- FindClusters(embryo.integrated, verbose = FALSE)
embryo.integrated <- RunUMAP(embryo.integrated, dims = 1:PCNum)

DimPlot(embryo.integrated, reduction = "umap", group.by = c("orig.ident"))
SpatialDimPlot(embryo.integrated)
Idents(embryo.integrated) <- "Embryo" 

cells_1 <- c("TCATAGCTTCCAGTAC-1_1",
             "GCGGTGGCCGACGCAT-1_1",
             "CTCAATGACCATTAAG-1_1",
             "GGCGTGTAATTACTTA-1_1",
             "TAAGTGGCCGCTCACC-1_1",
             "ATCCGGATGACGTAAC-1_1",
             "GGTCCTGCGAGTAGGT-1_1",
             "AAGTTCCTATTGCTGG-1_1",
             "GAACATTCGACACGTA-1_1",
             "CGGTATAGGATGAGCG-1_1",
             "AACCACTGCCATAGCC-1_1",
             "CAAGCCGCCATAGAGG-1_1"
)

Idents(embryo.integrated, cells_1) <- "Ovary"

cells_2 <- c(
  "AACCGCCAGACTACTT-1_2",
  "AATGGACATCCTACTC-1_2",
  "ATTAGTGCTTGGTCAT-1_2",
  "CGCAGATAGATGTTCT-1_2",
  "CGTAGTACGCATAATG-1_2",
  "CTACCAGACCTCTTAG-1_2",
  "GACAATTGTATGCTTC-1_2",
  "GCGGTGGCCGACGCAT-1_2",
  "GGCATTGCTAAGAATG-1_2",
  "GGTCTTGAACCGGCCA-1_2",
  "GGTTCACCGACTCACG-1_2",
  "TAACAGGTCCATACCA-1_2",
  "TAAGTGGCCGCTCACC-1_2",
  "TACTCCGAAGTAGAAT-1_2",
  "TAGTTCGAGTAGGAAG-1_2",
  "TCAGGAGCAACGGTCG-1_2",
  "TCATAGCTTCCAGTAC-1_2",
  "TCCGAGTACGACTTAG-1_2",
  "TCCGCCAGTGATTAAT-1_2",
  "TCGGCGCAGATAGATA-1_2",
  "TCTTACTCTGTATGCG-1_2"
)

Idents(embryo.integrated, cells_2) <- "Testis"
SpatialDimPlot(embryo.integrated)


embryo.integrated <- PrepSCTFindMarkers(embryo.integrated)
de_markers <- FindMarkers(embryo.integrated, ident.1 = "Ovary", ident.2 = "Testis")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Ovary_31_vs_Testis_29.csv")
