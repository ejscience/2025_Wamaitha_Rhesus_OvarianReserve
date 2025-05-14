library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)


#Load_19
emb_19 <- Load10X_Spatial('~/Desktop/New_Visium_Rm/emb019_outs/outs',filename = "filtered_feature_bc_matrix.h5", slice = "slice1")

#QC
VlnPlot(emb_19, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(emb_19, features = c("nCount_Spatial", "nFeature_Spatial"))
emb_19 <- emb_19[, emb_19$nFeature_Spatial > 100 & emb_19$nCount_Spatial > 100]
SpatialFeaturePlot(emb_19, features = c("nCount_Spatial", "nFeature_Spatial"))
dim(emb_19)

#SCT trasnform
emb_19 <- SCTransform(emb_19, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)
emb_19 <- RunPCA(emb_19, assay = "SCT", verbose = FALSE)
emb_19 <- FindNeighbors(emb_19, reduction = "pca", dims = 1:30)
emb_19 <- FindClusters(emb_19, verbose = FALSE)
emb_19 <- RunUMAP(emb_19, reduction = "pca", dims = 1:30)

p1 <- DimPlot(emb_19, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(emb_19, label = TRUE, label.size = 1, pt.size.factor = 1)
p1 + p2

cells_1 <- c(
  "AACTAGCCTTGCAATA-1",
  "AATAACGTAACTGTGA-1",
  "ACCAGTATGAAGTTCC-1",
  "ATCATTCACCAGATTG-1",
  "ATCCGTAGCTTCAGTA-1",
  "ATGTACTACCGACTTC-1",
  "ATTAGATAAGACGGCG-1",
  "CCGCGTCGTTGGTCTG-1",
  "CGCAACCATGCACCGA-1",
  "CTAATCTGGCCGTGTT-1",
  "GAGTTCTATCGTGCCT-1",
  "GCCTGGTACATGGTCG-1",
  "GCTCGGAGATGAATGT-1",
  "GTATGTTAGTTCTACC-1",
  "TAAGGACGCCTGGTTG-1",
  "TACGGATCAGCCATCC-1",
  "TAGGCATCTGCGGACC-1",
  "TCACCGTTGCTAATGG-1",
  "TCCGAACCTCAACCAT-1",
  "TCCGGCTCTATGGACA-1",
  "TCGTATCACGTTATGA-1",
  "TCTTGCACTATTCAGA-1",
  "TGATTAGAATGACATA-1",
  "TGGACGCCTGCGACCG-1",
  "TGGTTGCCTCCGGCTT-1",
  "GGCATAATCCTATCGC-1",
  "GAGCCTCACGTATGTA-1",
  "TACTCCTACATCAACT-1",
  "ATCGTACACCTATGCC-1",
  "GTCTTGACAATAGCAA-1",
  "TCCAGTACGTGCTATG-1",
  "TCAGCGTGACCTTACC-1",
  "CCGGCGAACCGTGCAC-1",
  "CTTGGTTAGAGTCATC-1"
)


Idents(emb_19) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_19, cells_1) <- "Testis"


SpatialDimPlot(emb_19, pt.size.factor = 1)

#Load_22
emb_22 <- Load10X_Spatial('~/Desktop/New_Visium_Rm/emb022_outs/outs/',filename = "filtered_feature_bc_matrix.h5", slice = "slice2")
#QC
dim(emb_22)
VlnPlot(emb_22, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(emb_22, features = c("nCount_Spatial", "nFeature_Spatial"))
emb_22 <- emb_22[, emb_22$nFeature_Spatial > 100 & emb_22$nCount_Spatial > 100]
SpatialFeaturePlot(emb_22, features = c("nCount_Spatial", "nFeature_Spatial"))
dim(emb_22)

#SCT trasnform
emb_22 <- SCTransform(emb_22, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)
emb_22 <- RunPCA(emb_22, assay = "SCT", verbose = FALSE)
emb_22 <- FindNeighbors(emb_22, reduction = "pca", dims = 1:30)
emb_22 <- FindClusters(emb_22, verbose = FALSE)
emb_22 <- RunUMAP(emb_22, reduction = "pca", dims = 1:30)

p1 <- DimPlot(emb_22, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(emb_22, label = TRUE, label.size = 3, pt.size.factor = 1)
p1 + p2

cells_6 <- c(
  "GCAGGCCATTCCTGCC-1",
  "GACAGTCACTTCGCTC-1",
  "CGACAATTATTCTGTC-1",
  "AAGCAGCGGACACGTG-1",
  "CTGTTCCAGAAGCTGA-1",
  "GCCAGTAGACGATACC-1",
  "ATCTGTACGCGTCTCC-1",
  "ATCTGTACGCGTCTCC-1",
  "TATCCACTTGCGTATC-1",
  "TCTTAGTCGTCTCGGC-1",
  "TATGACAATGCTTGAA-1",
  "TGCCTAAGTTAATTAT-1",
  "GTTCATCTTACATTCC-1",
  "TACGACGCGATTCAAG-1",
  "TGTCTTGCGACTAATT-1",
  "GTTCCAGTGCCTTACC-1",
  "TCAGCTAGGATGTACT-1",
  "GATTCCATACTGAACC-1",
  "ACCTGTGGTCACCTCG-1",
  "CTTACGTCATCCTTAT-1",
  "ATCCGGTATCACAAGT-1")


Idents(emb_22) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_22, cells_6) <- "Testis"


SpatialDimPlot(emb_22, pt.size.factor = 1)

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

emb_19@meta.data[["orig.ident"]] <- "emb_19"
emb_22@meta.data[["orig.ident"]] <- "emb_22"
emb_29@meta.data[["orig.ident"]] <- "emb_29"


#Integration
embryo_list <- list(Object_name1 = emb_19, 
                    Object_name2 = emb_22,
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
Idents(embryo.integrated) <- "ident" 


cells_1 <- c(
  "AACTAGCCTTGCAATA-1_1",
  "AATAACGTAACTGTGA-1_1",
  "ACCAGTATGAAGTTCC-1_1",
  "ATCATTCACCAGATTG-1_1",
  "ATCCGTAGCTTCAGTA-1_1",
  "ATGTACTACCGACTTC-1_1",
  "ATTAGATAAGACGGCG-1_1",
  "CCGCGTCGTTGGTCTG-1_1",
  "CGCAACCATGCACCGA-1_1",
  "CTAATCTGGCCGTGTT-1_1",
  "GAGTTCTATCGTGCCT-1_1",
  "GCCTGGTACATGGTCG-1_1",
  "GCTCGGAGATGAATGT-1_1",
  "GTATGTTAGTTCTACC-1_1",
  "TAAGGACGCCTGGTTG-1_1",
  "TACGGATCAGCCATCC-1_1",
  "TAGGCATCTGCGGACC-1_1",
  "TCACCGTTGCTAATGG-1_1",
  "TCCGAACCTCAACCAT-1_1",
  "TCCGGCTCTATGGACA-1_1",
  "TCGTATCACGTTATGA-1_1",
  "TCTTGCACTATTCAGA-1_1",
  "TGATTAGAATGACATA-1_1",
  "TGGACGCCTGCGACCG-1_1",
  "TGGTTGCCTCCGGCTT-1_1",
  "GGCATAATCCTATCGC-1_1",
  "GAGCCTCACGTATGTA-1_1",
  "TACTCCTACATCAACT-1_1",
  "ATCGTACACCTATGCC-1_1",
  "GTCTTGACAATAGCAA-1_1",
  "TCCAGTACGTGCTATG-1_1",
  "TCAGCGTGACCTTACC-1_1",
  "CCGGCGAACCGTGCAC-1_1",
  "CTTGGTTAGAGTCATC-1_1"
)


Idents(embryo.integrated, cells_1) <- "Testis_19"


cells_2 <- c(
  "GCAGGCCATTCCTGCC-1_2",
  "GACAGTCACTTCGCTC-1_2",
  "CGACAATTATTCTGTC-1_2",
  "AAGCAGCGGACACGTG-1_2",
  "CTGTTCCAGAAGCTGA-1_2",
  "GCCAGTAGACGATACC-1_2",
  "ATCTGTACGCGTCTCC-1_2",
  "ATCTGTACGCGTCTCC-1_2",
  "TATCCACTTGCGTATC-1_2",
  "TCTTAGTCGTCTCGGC-1_2",
  "TATGACAATGCTTGAA-1_2",
  "TGCCTAAGTTAATTAT-1_2",
  "GTTCATCTTACATTCC-1_2",
  "TACGACGCGATTCAAG-1_2",
  "TGTCTTGCGACTAATT-1_2",
  "GTTCCAGTGCCTTACC-1_2",
  "TCAGCTAGGATGTACT-1_2",
  "GATTCCATACTGAACC-1_2",
  "ACCTGTGGTCACCTCG-1_2",
  "CTTACGTCATCCTTAT-1_2",
  "ATCCGGTATCACAAGT-1_2"
)


Idents(embryo.integrated, cells_2) <- "Testis_22"


cells_3 <- c(
  "AACCGCCAGACTACTT-1_3",
  "AATGGACATCCTACTC-1_3",
  "ATTAGTGCTTGGTCAT-1_3",
  "CGCAGATAGATGTTCT-1_3",
  "CGTAGTACGCATAATG-1_3",
  "CTACCAGACCTCTTAG-1_3",
  "GACAATTGTATGCTTC-1_3",
  "GCGGTGGCCGACGCAT-1_3",
  "GGCATTGCTAAGAATG-1_3",
  "GGTCTTGAACCGGCCA-1_3",
  "GGTTCACCGACTCACG-1_3",
  "TAACAGGTCCATACCA-1_3",
  "TAAGTGGCCGCTCACC-1_3",
  "TACTCCGAAGTAGAAT-1_3",
  "TAGTTCGAGTAGGAAG-1_3",
  "TCAGGAGCAACGGTCG-1_3",
  "TCATAGCTTCCAGTAC-1_3",
  "TCCGAGTACGACTTAG-1_3",
  "TCCGCCAGTGATTAAT-1_3",
  "TCGGCGCAGATAGATA-1_3",
  "TCTTACTCTGTATGCG-1_3"
)

Idents(embryo.integrated, cells_3) <- "Testis_29"
SpatialDimPlot(embryo.integrated)

embryo.integrated <- PrepSCTFindMarkers(embryo.integrated)
de_markers <- FindMarkers(embryo.integrated, ident.1 = "Testis_19", ident.2 = "Testis_22")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Testis19_vs_Testis22.csv")

de_markers <- FindMarkers(embryo.integrated, ident.1 = "Testis_22", ident.2 = "Testis_29")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Testis22_vs_Testis29.csv")

de_markers <- FindMarkers(embryo.integrated, ident.1 = "Testis_19", ident.2 = "Testis_29")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Testis19_vs_Testis29.csv")

