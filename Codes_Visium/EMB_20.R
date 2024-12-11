library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

#Load
emb_20 <- Load10X_Spatial('~/Desktop/New_Visium_Rm/emb020_outs/outs',filename = "filtered_feature_bc_matrix.h5", slice = "slice1")

#QC
VlnPlot(emb_20, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(emb_20, features = c("nCount_Spatial", "nFeature_Spatial"))
emb_20 <- emb_20[, emb_20$nFeature_Spatial > 100 & emb_20$nCount_Spatial > 100]
SpatialFeaturePlot(emb_20, features = c("nCount_Spatial", "nFeature_Spatial"))
dim(emb_20)

#SCT trasnform
emb_20 <- SCTransform(emb_20, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)

emb_20 <- RunPCA(emb_20, assay = "SCT", verbose = FALSE)
emb_20 <- FindNeighbors(emb_20, reduction = "pca", dims = 1:30)
emb_20 <- FindClusters(emb_20, verbose = FALSE)
emb_20 <- RunUMAP(emb_20, reduction = "pca", dims = 1:30)

p1 <- DimPlot(emb_20, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(emb_20, label = TRUE, label.size = 3)
p1 + p2


#############mesonephoros##############
cells_1 <- c("AAGCAGCGGACACGTG-1",
             "AATAGTATAGTTGTCG-1",
             "AATTGCATTGAAGAAT-1",
             "ACGATGTTACCTCATT-1",
             "ACGTTAATACAATGAC-1",
             "ACTCAGTGAGCCTCCA-1",
             "AGAAGCAGCCTAAGCT-1",
             "AGACCGACCGCAGACA-1",
             "AGACTCATTGCGGACA-1",
             "AGCCGGCCTATTCCAG-1",
             "AGCGGTGAGTGACATG-1",
             "AGTCGCCACTCGGCGA-1",
             "AGTCGGTCTAATCTGT-1",
             "ATACATCAAGGATTGC-1",
             "ATAGTTCGCAGCCGTG-1",
             "ATCGCAATGCACTAAC-1",
             "ATCTGAGGAACGGTTG-1",
             "ATCTGTACGCGTCTCC-1",
             "ATGATATCTAGTCTTA-1",
             "ATTACATTGAAGTATG-1",
             "CACGTCCACGACCTAT-1",
             "CAGATTGCCTGTAAGT-1",
             "CATCAGATACTAATCG-1",
             "CCTAGCCTTCATCGGC-1",
             "CGACAATTATTCTGTC-1",
             "CGAGGAAGACTTAGAT-1",
             "CGAGTGCGGCTCGTGG-1",
             "CGCTCGGTACTGCGAA-1",
             "CGCTTATGGTTAGGAT-1",
             "CGGACTGTGAGTATAG-1",
             "CGTCAAGGCGGAATGT-1",
             "CGTCCATCGTCCAATA-1",
             "CTACATCTACTGGCTT-1",
             "CTCAATACGTTAACCG-1",
             "CTCACGAGTAGATGAG-1",
             "CTCCTACATTGTATCT-1",
             "CTGTTCCAGAAGCTGA-1",
             "CTTACGGTAAGCCGAA-1",
             "CTTAGCTACGCAAGGT-1",
             "CTTCGATACTCCGACC-1",
             "GAATAAGGAGAACTAA-1",
             "GAATGACACGACTAGG-1",
             "GACAGTCACTTCGCTC-1",
             "GAGGTGCATTCGGATC-1",
             "GATGCCTAGAGTCTGT-1",
             "GCAGGCCATTCCTGCC-1",
             "GCCAGTAGACGATACC-1",
             "GCCGGTATGCTCATGC-1",
             "GGAATCCTTAGACTGT-1",
             "GGACCGCTTCCTGGCA-1",
             "GGTAAGGCGACGATGA-1",
             "GTCTAAGATCGGACTA-1",
             "GTCTTGCGGCGATTCA-1",
             "GTGCAACGACGGTGTT-1",
             "GTGGTTGCCGCTTCAA-1",
             "GTTGACTAGAGCACGT-1",
             "TAATAGCGCTGTACTC-1",
             "TAGAGGCGATTGATTC-1",
             "TATCCACTTGCGTATC-1",
             "TCAAGGTATCAACTTA-1",
             "TCAGTCGCCAAGGTTG-1",
             "TCAGTTCTCCGGACCG-1",
             "TCCGTGTCTTACATAA-1",
             "TCGATCCAGTCTGTAG-1",
             "TCGTCCAGCATGGCTT-1",
             "TCTACTGTCTCAATGA-1",
             "TCTTAGTCGTCTCGGC-1",
             "TGATAGGAGCATCAAC-1",
             "TGATGGCCGGCATCGC-1",
             "TGCAGACAATACTGTT-1",
             "TGGCCTGGCACTGTGG-1",
             "TGGTCGCCGGCTTACT-1",
             "TGGTGCAGACCGAAGT-1")

cells_2 <- c("TATGAGCGATCCATTA-1", "GGATGGATCACTCATA-1", "CAGTGCCTTGCGGTCT-1", "AGGTGAGCAATGTCTT-1", "GGTAGCTTAGCCAATC-1",
             "TCATCTCCAACCACCT-1",
             "TATGTTGGTCCAGTGA-1",
             "TGTTGGCCAGACCTAC-1",
             "AGTTAAGACTGGCTCC-1",
             "CGACGCTCCGCACGTG-1", "GATGGACCTGGACGGA-1", "AGCGTGCTTCGCTCAT-1")

Idents(emb_20) <- "ident"  # Replace "ident" with the actual identity class in Seurat object
Idents(emb_20, cells_1) <- "Ovary"
Idents(emb_20, cells_2) <- "Mesonephoros"
SpatialDimPlot(emb_20, pt.size.factor = 1, crop = FALSE)

de_markers <- FindMarkers(emb_20, ident.1 = "Ovary")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Emb_20_Ovary_vs_ALL.csv")

de_markers <- FindMarkers(emb_20, ident.1 = "Ovary", ident.2 = "Mesonephoros")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Ovary_vs_Mesonephoros_emb_20.csv")
#####STOP####
