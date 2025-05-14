library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

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

de_markers <- FindMarkers(emb_22, ident.1 = "Testis")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_22_Testis_vs_ALL.csv")

de_markers <- FindMarkers(emb_22, ident.1 = "Testis", ident.2 = "Adrenal Gland")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_22_Testis_vs_Adrenal.csv")

de_markers <- FindMarkers(emb_22, ident.1 = "Testis_1", ident.2 = "Adrenal Gland")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_22_Testis_1_vs_Adrenal.csv")

de_markers <- FindMarkers(emb_22, ident.1 = "Testis", ident.2 = "Mesonephoros")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_22_Testis_ALL_vs_Mesonephoros_ALL.csv")

de_markers <- FindMarkers(emb_22, ident.1 = "Testis_1", ident.2 = "Mesonephoros_1")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_22_Testis_1_vs_Mesonephoros_1.csv")

de_markers <- FindMarkers(emb_22, ident.1 = "Testis_2", ident.2 = "Mesonephoros_2")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_22_Testis_2_vs_Mesonephoros_2.csv")

