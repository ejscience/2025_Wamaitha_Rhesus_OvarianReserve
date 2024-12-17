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

write.csv(emb_32@images[["slice1"]]@coordinates, file = "~/Desktop/Ovary_32_coordinates.txt")
write.csv(emb_32@active.ident, file = "~/Desktop/Ovary_32_idents.txt")

cells_1 <- c(
  "TAGGATGCACCGTTCA-1",
  "CACGATCATACAATTA-1",
  "TCAGTTCTCCGGACCG-1",
  "AGACCGACCGCAGACA-1",
  "CGAGGAAGACTTAGAT-1",
  "GCACCTTACTGGCTAC-1",
  "GCAGGCCATTCCTGCC-1",
  "AATAGTATAGTTGTCG-1",
  "CTGTCACTCGGCCTTA-1",
  "GAACAGATAACCTTAA-1",
  "GACAGTCACTTCGCTC-1",
  "GCCGGTATGCTCATGC-1",
  "AACTCAACCTTGACCA-1",
  "AAGCAGCGGACACGTG-1",
  "AGTGAGTATTCTAAGT-1",
  "CACGTCCACGACCTAT-1",
  "CTGTTCCAGAAGCTGA-1",
  "GTATCTGTCAACGGTG-1",
  "AGAAGCAGCCTAAGCT-1",
  "ATCTGTACGCGTCTCC-1",
  "ATGTGCGTCTGAACAG-1",
  "CGCCTTGTGAAGTCTA-1",
  "CTTCGATACTCCGACC-1",
  "GGCCTTCGACCAGGTA-1",
  "AACTTCGCTTAGTCAG-1",
  "CGCTCGGTACTGCGAA-1",
  "GGCTGACGCCGAATCA-1",
  "GTCTGATATCTGAATT-1",
  "TCTTCGCCATTCGAAT-1",
  "CAACTGCTGGAGTTAT-1",
  "GCGACACTAAGCTAAG-1",
  "TAGGTCTTGGCTCGTC-1",
  "AGTGCCAGATGCTCTG-1",
  "TAGTGACTCAATTGCT-1",
  "TGAGACAGAAGTGCTC-1"
)

cells_2 <- c(
  "GACAAGCTGAATCTGG-1",
  "ATTAAGAACAGCTCCG-1",
  "TGAAGTATGATCGTGC-1",
  "GAGTAGGATGGTAGTC-1",
  "CAAGACCACTATCGAA-1",
  "CTATGTAGGCCGAACG-1",
  "GTGACAAGATGCTCAT-1")

cells_3 <- c(
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

cells_4 <- c(
  "AGCCTCCGAACTACCA-1",
  "GGACAGTGCGGCAACA-1",
  "ACCAGTATGAAGTTCC-1",
  "CGCAACCATGCACCGA-1",
  "GCCACTATAGACCATA-1",
  "GCTCGGAGATGAATGT-1",
  "GGCATAATCCTATCGC-1",
  "TAGACCGGCCAGTACG-1",
  "TCCGGCTCTATGGACA-1",
  "TCGTATCACGTTATGA-1",
  "AACTAGCACTGTTCGA-1",
  "AGCCTCCGAACTACCA-1",
  "ATGTACTACCGACTTC-1",
  "CCGCGTCGTTGGTCTG-1",
  "CTGGAATGGATCAGCA-1",
  "GAGCCTCACGTATGTA-1",
  "TCCAGTACGTGCTATG-1",
  "TGCCAATCCATGTCTA-1",
  "AACTAGCCTTGCAATA-1",
  "TCTGGATATTAACCTA-1",
  "CTAATCTGGCCGTGTT-1",
  "TAGGCATCTGCGGACC-1")

cells_5 <- c("ATTAGATAAGACGGCG-1",
"GCCTGGTACATGGTCG-1",
"GTATGTTAGTTCTACC-1",
"TGATTAGAATGACATA-1",
"CTTGGTTAGAGTCATC-1",
"TCACCGTTGCTAATGG-1",
"TCCGAACCTCAACCAT-1",
"AATAACGTAACTGTGA-1",
"GAGTTCTATCGTGCCT-1",
"ATCCGTAGCTTCAGTA-1",
"ATCATTCACCAGATTG-1")

cells_6 <- c("ACCAAGTGGCTGTATG-1",
             "ATCGTACACCTATGCC-1",
             "GTTAACTCCGACCACT-1",
             "TAAGGACGCCTGGTTG-1",
             "CGAGGAGTCTAAGGAT-1",
             "CTGGTGCGTTCACGAG-1",
             "GCGAGCGGACCATAAC-1",
             "GTCTTGACAATAGCAA-1",
             "TACTCCTACATCAACT-1",
             "TCAGCGTGACCTTACC-1")

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

cells_8 <- c(
  "TAGGATGCACCGTTCA-1",
  "CACGATCATACAATTA-1",
  "TCAGTTCTCCGGACCG-1",
  "AGACCGACCGCAGACA-1",
  "CGAGGAAGACTTAGAT-1",
  "GCACCTTACTGGCTAC-1",
  "GCAGGCCATTCCTGCC-1",
  "AATAGTATAGTTGTCG-1",
  "CTGTCACTCGGCCTTA-1",
  "GAACAGATAACCTTAA-1",
  "GACAGTCACTTCGCTC-1",
  "GCCGGTATGCTCATGC-1",
  "AACTCAACCTTGACCA-1",
  "AAGCAGCGGACACGTG-1",
  "AGTGAGTATTCTAAGT-1",
  "CACGTCCACGACCTAT-1",
  "CTGTTCCAGAAGCTGA-1",
  "GTATCTGTCAACGGTG-1",
  "AGAAGCAGCCTAAGCT-1",
  "ATCTGTACGCGTCTCC-1",
  "ATGTGCGTCTGAACAG-1",
  "CGCCTTGTGAAGTCTA-1",
  "CTTCGATACTCCGACC-1",
  "GGCCTTCGACCAGGTA-1",
  "AACTTCGCTTAGTCAG-1",
  "CGCTCGGTACTGCGAA-1",
  "GGCTGACGCCGAATCA-1",
  "GTCTGATATCTGAATT-1",
  "TCTTCGCCATTCGAAT-1",
  "CAACTGCTGGAGTTAT-1",
  "GCGACACTAAGCTAAG-1",
  "TAGGTCTTGGCTCGTC-1",
  "AGTGCCAGATGCTCTG-1",
  "TAGTGACTCAATTGCT-1",
  "TGAGACAGAAGTGCTC-1",
  "AGCCTCCGAACTACCA-1",
  "GGACAGTGCGGCAACA-1",
  "ACCAGTATGAAGTTCC-1",
  "CGCAACCATGCACCGA-1",
  "GCCACTATAGACCATA-1",
  "GCTCGGAGATGAATGT-1",
  "GGCATAATCCTATCGC-1",
  "TAGACCGGCCAGTACG-1",
  "TCCGGCTCTATGGACA-1",
  "TCGTATCACGTTATGA-1",
  "AACTAGCACTGTTCGA-1",
  "AGCCTCCGAACTACCA-1",
  "ATGTACTACCGACTTC-1",
  "CCGCGTCGTTGGTCTG-1",
  "CTGGAATGGATCAGCA-1",
  "GAGCCTCACGTATGTA-1",
  "TCCAGTACGTGCTATG-1",
  "TGCCAATCCATGTCTA-1",
  "AACTAGCCTTGCAATA-1",
  "TCTGGATATTAACCTA-1",
  "CTAATCTGGCCGTGTT-1",
  "TAGGCATCTGCGGACC-1"
)

cells_9 <- c("ACCAAGTGGCTGTATG-1",
             "ATCGTACACCTATGCC-1",
             "GTTAACTCCGACCACT-1",
             "TAAGGACGCCTGGTTG-1",
             "CGAGGAGTCTAAGGAT-1",
             "CTGGTGCGTTCACGAG-1",
             "GCGAGCGGACCATAAC-1",
             "GTCTTGACAATAGCAA-1",
             "TACTCCTACATCAACT-1",
             "TCAGCGTGACCTTACC-1",
             "GACAAGCTGAATCTGG-1",
             "ATTAAGAACAGCTCCG-1",
             "TGAAGTATGATCGTGC-1",
             "GAGTAGGATGGTAGTC-1",
             "CAAGACCACTATCGAA-1",
             "CTATGTAGGCCGAACG-1",
             "GTGACAAGATGCTCAT-1")

Idents(emb_32) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_32, cells_1) <- "Mesonephoros_1"
Idents(emb_32, cells_2) <- "Adrenal_1"
Idents(emb_32, cells_3) <- "Ovary_1"
Idents(emb_32, cells_4) <- "Mesonephoros_2"
Idents(emb_32, cells_5) <- "Ovary_2"
Idents(emb_32, cells_6) <- "Adrenal_2"
Idents(emb_32, cells_7) <- "Ovaries"
Idents(emb_32, cells_8) <- "Mesonephoros"
Idents(emb_32, cells_9) <- "Adrenal Gland"

SpatialDimPlot(emb_32, pt.size.factor = 1)

de_markers <- FindMarkers(emb_32, ident.1 = "Ovaries")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_32_Ovaries_vs_ALL.csv")


de_markers <- FindMarkers(emb_32, ident.1 = "Ovaries", ident.2 = "Mesonephoros")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_32_Ovary_ALL_vs_Mesonephoros_ALL.csv")

de_markers <- FindMarkers(emb_32, ident.1 = "Ovaries", ident.2 = "Adrenal Gland")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_32_Ovary_ALL_vs_Adrenal_ALL.csv")

de_markers <- FindMarkers(emb_32, ident.1 = "Ovary_1", ident.2 = "Mesonephoros_1")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_32_Ovary_1_vs_Mesonephoros_1.csv")

de_markers <- FindMarkers(emb_32, ident.1 = "Ovary_2", ident.2 = "Mesonephoros_2")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_32_Ovary_2_vs_Mesonephoros_2.csv")

de_markers <- FindMarkers(emb_32, ident.1 = "Ovary_1", ident.2 = "Adrenal_1")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_32_Ovary_1_vs_Adrenal_1.csv")

de_markers <- FindMarkers(emb_32, ident.1 = "Ovary_2", ident.2 = "Adrenal_2")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_32_Ovary_2_vs_Adrenal_2.csv")



