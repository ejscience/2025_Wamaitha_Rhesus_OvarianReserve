library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

#Load
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

write.csv(x@images[["slice1"]]@coordinates, file = "~/Desktop/Ovary_19_coordinates.txt")

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

cells_2 <- c("CGGCAACCACTAGCAT-1",
             "AAGTGGCGAGTGTACA-1",
             "GTACTTAACTTCAGCG-1",
             "GCAACATACATCAGCA-1",
             "AAGGAGGCACTCAATT-1",
             "GGTACCAGGTAGTGTT-1",
             "AACGAAGCGTGGAAGT-1",
             "TCAGCTACCAACTCCG-1",
             "GCTTGCCATAGTGGAT-1",
             "GCGAGCGGACCATAAC-1",
             "CTTCCGTCGCTGCATC-1",
             "TAGACCGGCCAGTACG-1",
             "AGCCTCCGAACTACCA-1",
             "CAGGACATGGTAATGC-1",
             "ACGGCTTCCATCGGAA-1",
             "GAACTAACTCTACGTG-1",
             "CACGCCACAATCGCAC-1",
             "AGTACCGTTGTACGCC-1",
             "GTGAACATTCAACGAT-1",
             "ATTCGCTCGCGGAGGC-1",
             "CCGTGAAGGATAACAG-1",
             "CACGTAGGTAGGCTGA-1",
             "CAGCTCCTCGGATTCT-1",
             "GGTAGCCATCAATCCG-1",
             "GCTTCTTGGACGAGCG-1",
             "CACGTCAAGTGCGATC-1",
             "TAACGAGTTGCCAACC-1",
             "AGCGTGCGCCTGATCG-1",
             "ACGGCTTGTTCAATCG-1")



Idents(emb_19) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_19, cells_1) <- "Testis"
Idents(emb_19, cells_2) <- "Mesonephoros"

SpatialDimPlot(emb_19, pt.size.factor = 1)

de_markers <- FindMarkers(emb_19, ident.1 = "Testis")
write.csv(de_markers, file = "~/Desktop//NEW_Visium_Rm/Results/Emb_19_Testis_vs_ALL.csv")

de_markers <- FindMarkers(emb_19, ident.1 = "Testis", ident.2 = "Mesonephoros")
write.csv(de_markers, file = "~/Desktop//NEW_Visium_Rm/Results/emb_19_Testis_vs_Mesonephoros.csv")
