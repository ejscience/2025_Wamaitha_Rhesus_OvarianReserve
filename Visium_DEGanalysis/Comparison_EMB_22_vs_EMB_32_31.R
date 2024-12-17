library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

#Load
emb_32 <- Load10X_Spatial('~/Desktop/New_Visium_Rm/emb032_outs/outs/',filename = "filtered_feature_bc_matrix.h5", slice = "slice1")

#QC
VlnPlot(emb_32, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(emb_32, features = c("nCount_Spatial", "nFeature_Spatial"))
emb_32 <- emb_32[, emb_32$nFeature_Spatial > 100 & emb_32$nCount_Spatial > 100]
SpatialFeaturePlot(emb_32, features = c("nCount_Spatial", "nFeature_Spatial"))
dim(emb_32)

#SCT trasnform
emb_32 <- SCTransform(emb_32, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)
emb_32 <- RunPCA(emb_32, assay = "SCT", verbose = FALSE)
emb_32 <- FindNeighbors(emb_32, reduction = "pca", dims = 1:30)
emb_32 <- FindClusters(emb_32, verbose = FALSE)
emb_32 <- RunUMAP(emb_32, reduction = "pca", dims = 1:30)

p1 <- DimPlot(emb_32, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(emb_32, label = TRUE, label.size = 1, pt.size.factor = 1)
p1 + p2

cells_7 <- c("ATTAGATAAGACGGCG-1",
             "GCCTGGTACATGGTCG-1",
             "GTATGTTAGTTCTACC-1",
             "TGATTAGAATGACATA-1",
             "CTTGGTTAGAGTCATC-1",
             "TCACCGTTGCTAATGG-1",
             "TCCGAACCTCAACCAT-1",
             "AATAACGTAACTGTGA-1",
             "GAGTTCTATCGTGCCT-1",
             "ATCCGTAGCTTCAGTA-1",
             "ATCATTCACCAGATTG-1",
             "CTTAGCTACGCAAGGT-1",
             "TATGACAATGCTTGAA-1",
             "TGCCTAAGTTAATTAT-1",
             "ATACGCTGAGATCAGA-1",
             "TATCCACTTGCGTATC-1",
             "TCTTAGTCGTCTCGGC-1",
             "CGACAATTATTCTGTC-1",
             "CGTCAAGGCGGAATGT-1",
             "GAGGTGCATTCGGATC-1",
             "GCCAGTAGACGATACC-1")


Idents(emb_32) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_32, cells_7) <- "Ovaries"


SpatialDimPlot(emb_32, pt.size.factor = 1)


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

#Load
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
p2 <- SpatialDimPlot(emb_22, label = TRUE, label.size = 3, pt.size.factor = 0.1)
p1 + p2

LinkedDimPlot(emb_22, alpha = c(0.001, 0.001), combine = FALSE)
write.csv(emb_22@images[["slice2"]]@coordinates, file = "~/Desktop/Ovary_22_coordinates.txt")
write.csv(emb_31@active.ident, file = "~/Desktop/Ovary_31_idents.txt")

#GTTCATCTTACATTCC-1 T2
cells_1 <- c(
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
  "TGCCTAAGTTAATTAT-1")

cells_2 <- c(  
  "GTTCATCTTACATTCC-1",
  "TACGACGCGATTCAAG-1",
  "TGTCTTGCGACTAATT-1",
  "GTTCCAGTGCCTTACC-1",
  "TCAGCTAGGATGTACT-1",
  "GATTCCATACTGAACC-1",
  "ACCTGTGGTCACCTCG-1",
  "CTTACGTCATCCTTAT-1",
  "ATCCGGTATCACAAGT-1"
)

#"GATGAGTGAACGGCTG-1" ME1
cells_3 <- c("TAATTACCGGTGCGAG-1",
             "CGGTCCATTACAACTA-1",
             "GTCAGCCTCCGAACCT-1",
             "AACCGCCAGACTACTT-1",
             "TACTCCGAAGTAGAAT-1",
             "TGAACGGACACGTATC-1",
             "TCCAACTGAATAAGGT-1",
             "CGATAAGAATTCCAGT-1",
             "TCGGCGCAGATAGATA-1",
             "AGACTCATGACCTGGC-1",
             "GCCTTCTGTCACGCTC-1",
             "TAGACCGCGGAAGGAC-1",
             "ACTCTACTTCCGCATT-1",
             "GTGGCGCACCTGTTCC-1",
             "GTCAGTAGGTGACCGG-1",
             "CACTACTGCGGAATTG-1",
             "ATAGACATGGCTATGG-1",
             "TCTGAGTCAGATCTTG-1",
             "AGCACGCTCGATCTAA-1",
             "CTTAGTGGTACACGAC-1",
             "GTCACTTATCCAAGAC-1",
             "ACGGTTCGTACAGTTG-1",
             "AAGACTGTCATAGTGC-1",
             "CCACGTGCCACCTCAA-1")

cells_4 <- c(
  "GATGAGTGAACGGCTG-1",
  "CTGCCGTGAAGTCGCA-1",
  "CTTCGATAACATTGGT-1",
  "GCGAATCGAGAACACG-1",
  "CATCAGATACTAATCG-1",
  "TGCATTGCTGTCGGCG-1",
  "TAGGATGCACCGTTCA-1",
  "TGCAACGTTACCTGGT-1",
  "CGCTTATGGTTAGGAT-1",
  "TCAGTTCTCCGGACCG-1",
  "CACGATCATACAATTA-1",
  "CAGATTCCGTATGCGA-1",
  "AGACCGACCGCAGACA-1",
  "CGAGGAAGACTTAGAT-1",
  "GCACCTTACTGGCTAC-1")

cells_5 <- c(
  "GATGAGTGAACGGCTG-1",
  "CTGCCGTGAAGTCGCA-1",
  "CTTCGATAACATTGGT-1",
  "GCGAATCGAGAACACG-1",
  "CATCAGATACTAATCG-1",
  "TGCATTGCTGTCGGCG-1",
  "TAGGATGCACCGTTCA-1",
  "TGCAACGTTACCTGGT-1",
  "CGCTTATGGTTAGGAT-1",
  "TCAGTTCTCCGGACCG-1",
  "CACGATCATACAATTA-1",
  "CAGATTCCGTATGCGA-1",
  "AGACCGACCGCAGACA-1",
  "CGAGGAAGACTTAGAT-1",
  "GCACCTTACTGGCTAC-1",
  "TAATTACCGGTGCGAG-1",
  "CGGTCCATTACAACTA-1",
  "GTCAGCCTCCGAACCT-1",
  "AACCGCCAGACTACTT-1",
  "TACTCCGAAGTAGAAT-1",
  "TGAACGGACACGTATC-1",
  "TCCAACTGAATAAGGT-1",
  "CGATAAGAATTCCAGT-1",
  "TCGGCGCAGATAGATA-1",
  "AGACTCATGACCTGGC-1",
  "GCCTTCTGTCACGCTC-1",
  "TAGACCGCGGAAGGAC-1",
  "ACTCTACTTCCGCATT-1",
  "GTGGCGCACCTGTTCC-1",
  "GTCAGTAGGTGACCGG-1",
  "CACTACTGCGGAATTG-1",
  "ATAGACATGGCTATGG-1",
  "TCTGAGTCAGATCTTG-1",
  "AGCACGCTCGATCTAA-1",
  "CTTAGTGGTACACGAC-1",
  "GTCACTTATCCAAGAC-1",
  "ACGGTTCGTACAGTTG-1",
  "AAGACTGTCATAGTGC-1",
  "CCACGTGCCACCTCAA-1")

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

cells_7 <- c(
  "AGACTCATTGCGGACA-1",
  "TGGTCGCCGGCTTACT-1",
  "CAGATTGCCTGTAAGT-1",
  "GTGGTTGCCGCTTCAA-1",
  "GGTAAGGCGACGATGA-1",
  "TAGAGGCGATTGATTC-1",
  "AGCGGTGAGTGACATG-1",
  "TCGTCCAGCATGGCTT-1",
  "CGTCAAGGCGGAATGT-1",
  "CGAGTGCGGCTCGTGG-1",
  "GATGCCTAGAGTCTGT-1")

Idents(emb_22) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_22, cells_1) <- "Testis_1"
Idents(emb_22, cells_2) <- "Testis_2"
Idents(emb_22, cells_3) <- "Mesonephoros_1"
Idents(emb_22, cells_4) <- "Mesonephoros_2"
Idents(emb_22, cells_6) <- "Testis"
Idents(emb_22, cells_5) <- "Mesonephoros"
Idents(emb_22, cells_7) <- "Adrenal Gland"

SpatialDimPlot(emb_22, pt.size.factor = 1)


emb_22@meta.data[["orig.ident"]] <- "emb_22"
emb_31@meta.data[["orig.ident"]] <- "emb_31"
emb_32@meta.data[["orig.ident"]] <- "emb_32"

#Integration
embryo_list <- list(Object_name1 = emb_31, 
                    Object_name2 = emb_32,
                    Object_name3 = emb_22)

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


cells_1 <- c("ATTAGATAAGACGGCG-1_2",
             "GCCTGGTACATGGTCG-1_2",
             "GTATGTTAGTTCTACC-1_2",
             "TGATTAGAATGACATA-1_2",
             "CTTGGTTAGAGTCATC-1_2",
             "TCACCGTTGCTAATGG-1_2",
             "TCCGAACCTCAACCAT-1_2",
             "AATAACGTAACTGTGA-1_2",
             "GAGTTCTATCGTGCCT-1_2",
             "ATCCGTAGCTTCAGTA-1_2",
             "ATCATTCACCAGATTG-1_2",
             "CTTAGCTACGCAAGGT-1_2",
             "TATGACAATGCTTGAA-1_2",
             "TGCCTAAGTTAATTAT-1_2",
             "ATACGCTGAGATCAGA-1_2",
             "TATCCACTTGCGTATC-1_2",
             "TCTTAGTCGTCTCGGC-1_2",
             "CGACAATTATTCTGTC-1_2",
             "CGTCAAGGCGGAATGT-1_2",
             "GAGGTGCATTCGGATC-1_2",
             "GCCAGTAGACGATACC-1_2")

Idents(embryo.integrated, cells_1) <- "Gonad_2"


cells_2 <- c("TCATAGCTTCCAGTAC-1_1",
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

Idents(embryo.integrated, cells_2) <- "Gonad_1"

cells_3 <- c(
  "GCAGGCCATTCCTGCC-1_3",
  "GACAGTCACTTCGCTC-1_3",
  "CGACAATTATTCTGTC-1_3",
  "AAGCAGCGGACACGTG-1_3",
  "CTGTTCCAGAAGCTGA-1_3",
  "GCCAGTAGACGATACC-1_3",
  "ATCTGTACGCGTCTCC-1_3",
  "ATCTGTACGCGTCTCC-1_3",
  "TATCCACTTGCGTATC-1_3",
  "TCTTAGTCGTCTCGGC-1_3",
  "TATGACAATGCTTGAA-1_3",
  "TGCCTAAGTTAATTAT-1_3",
  "GTTCATCTTACATTCC-1_3",
  "TACGACGCGATTCAAG-1_3",
  "TGTCTTGCGACTAATT-1_3",
  "GTTCCAGTGCCTTACC-1_3",
  "TCAGCTAGGATGTACT-1_3",
  "GATTCCATACTGAACC-1_3",
  "ACCTGTGGTCACCTCG-1_3",
  "CTTACGTCATCCTTAT-1_3",
  "ATCCGGTATCACAAGT-1_3")


Idents(embryo.integrated, cells_3) <- "Gonad_3"

SpatialDimPlot(embryo.integrated)

embryo.integrated <- PrepSCTFindMarkers(embryo.integrated)
de_markers <- FindMarkers(embryo.integrated, ident.1 = "Gonad_3", ident.2 = "Gonad_2")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Gonad_22_vs_Gonad_32.csv")

de_markers <- FindMarkers(embryo.integrated, ident.1 = "Gonad_3", ident.2 = "Gonad_1")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Gonad_22_vs_Gonad_31.csv")
