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

write.csv(emb_31@images[["slice1"]]@coordinates, file = "~/Desktop/Ovary_31_coordinates.txt")
write.csv(emb_31@active.ident, file = "~/Desktop/Ovary_31_idents.txt")


cells_1 <- c(
  "ATCCGGATGACGTAAC-1",
  "GGTCCTGCGAGTAGGT-1",
  "AAGTTCCTATTGCTGG-1",
  "GAACATTCGACACGTA-1",
  "CGGTATAGGATGAGCG-1",
  "AACCACTGCCATAGCC-1",
  "CAAGCCGCCATAGAGG-1"
)

cells_2 <- c("ACCGGCACCATACTGC-1",
            "ATCCTCGTTCCAGGCC-1",
            "CTTCCAATTGAAGATT-1",
            "TCACATCCTTCTCAGT-1",
            "TAGATGGCTGTATGCC-1",
            "GTTGGTCCGCCAATCT-1",
            "ATAGATTGAGTTCGAG-1",
            "GGAACTTCGACATGTC-1",
            "CCGCGTTACCTTGGCG-1",
            "AGTTGTCCACAGACCA-1",
            "TATGAGGATCTACTGG-1",
            "CTCGTCGCGATTAATT-1",
            "CATGCTCATCTGAACG-1",
            "TGTCAGAAGTTGTAAC-1")

cells_3 <- c("ACGACTGAAGACTTAC-1",
             "GTCAGTTCTGCCTGGA-1",
             "TCATTGTTCAATGCGT-1",
             "TAACGCGTAATATGGC-1"
             )

cells_4 <- c("ACAGGAACTACAGAAC-1",
             "GCGGCAAGTGTTATGT-1",
             "GTTCATCTTACATTCC-1")

cells_5 <- c("CGATGTCCGACTACAG-1",
             "CGTCGGCAGACGGTAT-1",
             "GGCATTGCTAAGAATG-1",
             "TCAGGAGCAACGGTCG-1",
             "GTCTCATGACTAGTTG-1",
             "TCCGCCAGTGATTAAT-1",
             "AATGGACATCCTACTC-1",
             "AGACGGCCTTACGGAG-1",
             "CGTAGTACGCATAATG-1",
             "GGTCTTGAACCGGCCA-1",
             "TCTTACTCTGTATGCG-1")

cells_6 <- c("TCATAGCTTCCAGTAC-1",
             "GCGGTGGCCGACGCAT-1",
             "CTCAATGACCATTAAG-1",
             "GGCGTGTAATTACTTA-1",
             "TAAGTGGCCGCTCACC-1"
             )

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

cells_8 <- c("ACCGGCACCATACTGC-1",
             "ATCCTCGTTCCAGGCC-1",
             "CTTCCAATTGAAGATT-1",
             "TCACATCCTTCTCAGT-1",
             "TAGATGGCTGTATGCC-1",
             "GTTGGTCCGCCAATCT-1",
             "ATAGATTGAGTTCGAG-1",
             "GGAACTTCGACATGTC-1",
             "CCGCGTTACCTTGGCG-1",
             "AGTTGTCCACAGACCA-1",
             "TATGAGGATCTACTGG-1",
             "CTCGTCGCGATTAATT-1",
             "CATGCTCATCTGAACG-1",
             "TGTCAGAAGTTGTAAC-1",
             "CGATGTCCGACTACAG-1",
             "CGTCGGCAGACGGTAT-1",
             "GGCATTGCTAAGAATG-1",
             "TCAGGAGCAACGGTCG-1",
             "GTCTCATGACTAGTTG-1",
             "TCCGCCAGTGATTAAT-1",
             "AATGGACATCCTACTC-1",
             "AGACGGCCTTACGGAG-1",
             "CGTAGTACGCATAATG-1",
             "GGTCTTGAACCGGCCA-1",
             "TCTTACTCTGTATGCG-1")

cells_9 <- c("ACGACTGAAGACTTAC-1",
             "GTCAGTTCTGCCTGGA-1",
             "TCATTGTTCAATGCGT-1",
             "TAACGCGTAATATGGC-1","ACAGGAACTACAGAAC-1",
             "GCGGCAAGTGTTATGT-1",
             "GTTCATCTTACATTCC-1"
)


Idents(emb_31) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_31, cells_1) <- "Ovary_1"
Idents(emb_31, cells_2) <- "Mesonephoros_1"
Idents(emb_31, cells_3) <- "Adrenal_1"
Idents(emb_31, cells_4) <- "Adrenal_2"
Idents(emb_31, cells_5) <- "Mesonephoros_2"
Idents(emb_31, cells_6) <- "Ovary_2"
Idents(emb_31, cells_7) <- "Ovaries"
Idents(emb_31, cells_8) <- "Mesonephoros"
Idents(emb_31, cells_9) <- "Adrenal Gland"

SpatialDimPlot(emb_31, pt.size.factor = 1)

de_markers <- FindMarkers(emb_31, ident.1 = "Ovaries")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Emb_31_Ovaryies_vs_ALL.csv")

de_markers <- FindMarkers(emb_31, ident.1 = "Ovaries", ident.2 = "Mesonephoros")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/emb_31_Ovary_ALL_vs_Mesonephoros_ALL.csv")

de_markers <- FindMarkers(emb_31, ident.1 = "Ovaries", ident.2 = "Adrenal Gland")
write.csv(de_markers, file = "~/Desktop/emb_31_Ovary_ALL_vs_Adrenal_ALL.csv")

de_markers <- FindMarkers(emb_31, ident.1 = "Ovary_1", ident.2 = "Mesonephoros_1")
write.csv(de_markers, file = "~/Desktop/emb_31_Ovary_1_vs_Mesonephoros_1.csv")

de_markers <- FindMarkers(emb_31, ident.1 = "Ovary_2", ident.2 = "Mesonephoros_2")
write.csv(de_markers, file = "~/Desktop/emb_31_Ovary_2_vs_Mesonephoros_2.csv")

de_markers <- FindMarkers(emb_31, ident.1 = "Ovary_1", ident.2 = "Adrenal_1")
write.csv(de_markers, file = "~/Desktop/emb_31_Ovary_1_vs_Adrenal_1.csv")

de_markers <- FindMarkers(emb_31, ident.1 = "Ovary_2", ident.2 = "Adrenal_2")
write.csv(de_markers, file = "~/Desktop/emb_31_Ovary_2_vs_Adrenal_2.csv")

