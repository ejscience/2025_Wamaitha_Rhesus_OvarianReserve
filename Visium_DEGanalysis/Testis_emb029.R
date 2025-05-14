library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)


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


LinkedDimPlot(emb_29)

x <- subset(emb_29, idents = c(3, 7, 13))
SpatialDimPlot(x, label = TRUE, label.size = 3, pt.size.factor = 2)
write.csv(x@images[["slice3"]]@coordinates, file = "~/Desktop/Ovary_19_coordinates.txt")

LinkedDimPlot(x)


cells_1 <- c(
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
  "ATTAGTGCTTGGTCAT-1"
  )
  
cells_2 <- c(
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

cells_3 <- c(
             "AAGCCAGGCCGACGTG-1",
             "AAGTCTTAAGTTGCCA-1",
             "ACCAAGTGGCTGTATG-1",
             "AGCGAGAACAGATTCT-1",
             "AGTTGTACTTAATGGA-1",
             "ATCCGTAGCTTCAGTA-1",
             "CCGAGGATGACTCGTC-1",
             "CCGGCGAACCGTGCAC-1",
             "CGATGCACCAGGCCTC-1",
             "CGGTTATCTGTCTTAA-1",
             "CTATGTGCGGTTGACT-1",
             "CTGGAATAGTGCAAGT-1",
             "CTTAGTCTATCGAGAC-1",
             "GACAACGGACTCAGCT-1",
             "GACCACATTCTTGAAT-1",
             "GATGTAGGAGTTGGCC-1",
             "GCAAGATAGGTTAACT-1",
             "GCTCATACAATTAGAC-1",
             "GCCTTGGCGAACAATC-1",
             "TAACGTCACTATAATC-1",
             "GTACCTGAAGTGAATG-1",
             "TGGTAGCATCAAGTCT-1",
             "TGAGGTAGGACTACTG-1",
             "TCAGCGTGACCTTACC-1",
             "GTATGTTAGTTCTACC-1",
             "CTTGGTTAGAGTCATC-1",
             "TGCAGTTATTCCACTC-1",
             "GCCTTGCCACTTCTGA-1",
             "TCACAAGCTCGAGGCA-1",
             "GTGAAGATGAGTGCCT-1",
             "CTGCTCATTAGAGTGC-1",
             "TCTAGGCGGAAGGCAC-1",
             "GGCGGAAGTCTTGTAG-1",
             "AACTGCCTCGATAGGT-1",
             "CTCAATGACCATTAAG-1",
             "TGGCCTGACATACTCA-1",
             "GGCGTGTAATTACTTA-1")

#"GATGAGTGAACGGCTG-1" ME1
cells_4 <- c("AACGAATTGACCGGTT-1",
             "GGAGGTCGTAGGAGCT-1",
             "GTACAGATGACATCGC-1",
             "GTCAGCCTCCGAACCT-1",
             "GTCATCACTATTGTGC-1",
             "GTTGAGCTTCGTCCGG-1",
             "TAATTACCGGTGCGAG-1",
             "TAGGACATCAGCATGG-1",
             "TAGTGAATGTTCAGTG-1",
             "TAGTGCAGTGCGGTAG-1",
             "TAGTTCCGATAGAATC-1",
             "TCCGAATCATGTTGTC-1",
             "TGGATCATTGTAACGG-1",
             "CAGAGGCACTTAAGTC-1",
             "TGTAATCATATGAGAC-1",
             "TGTTCAGAACGGTGTA-1",
             "TGTTCTGCTCTGTCGT-1",
             "GGTCTTCCGCACATCC-1",
             "GATGAGTGAACGGCTG-1",
             "CTTAGTGGTACACGAC-1",
             "TCTGAGTCAGATCTTG-1",
             "ACTCTACTTCCGCATT-1",
             "AGACTCATGACCTGGC-1",
             "TAGACCGCGGAAGGAC-1",
             "AGCACGCTCGATCTAA-1",
             "CGCCATGAGGCGATGA-1",
             "TGAACGGATCTGGTAC-1",
             "CGAGTCATCCGCACCA-1",
             "CGGTCCATTACAACTA-1",
             "TGTCTTGGAATTGGCA-1",
             "CGTTGGTCGTCACTGT-1",
             "CGGAGCGGTAAGACGC-1",
             "CGTACATCACTGGATG-1",
             "GCCTTCTGTCACGCTC-1",
             "GCACCTAGCCGAATCT-1",
             "CAATAACTACATTAAC-1",
             "CACCATAAGCTGGCTG-1",
             "CTCCACGCTAGCTAGA-1",
             "GCATGCTGAGACAGGT-1",
             "GAGGTAACTGCGAGAC-1",
             "ACGCTGCGATTAACGA-1",
             "ATCTGGATCCTCAGAC-1",
             "CCTACACCGACATTAC-1",
             "CGCTCGCAATTGGCGA-1"
             )

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

cells_6 <- c(
  "AAGCCAGGCCGACGTG-1",
  "AAGTCTTAAGTTGCCA-1",
  "ACCAAGTGGCTGTATG-1",
  "AGCGAGAACAGATTCT-1",
  "AGTTGTACTTAATGGA-1",
  "ATCCGTAGCTTCAGTA-1",
  "CCGAGGATGACTCGTC-1",
  "CCGGCGAACCGTGCAC-1",
  "CGATGCACCAGGCCTC-1",
  "CGGTTATCTGTCTTAA-1",
  "CTATGTGCGGTTGACT-1",
  "CTGGAATAGTGCAAGT-1",
  "CTTAGTCTATCGAGAC-1",
  "GACAACGGACTCAGCT-1",
  "GACCACATTCTTGAAT-1",
  "GATGTAGGAGTTGGCC-1",
  "GCAAGATAGGTTAACT-1",
  "GCTCATACAATTAGAC-1",
  "GCCTTGGCGAACAATC-1",
  "TAACGTCACTATAATC-1",
  "GTACCTGAAGTGAATG-1",
  "TGGTAGCATCAAGTCT-1",
  "TGAGGTAGGACTACTG-1",
  "TCAGCGTGACCTTACC-1",
  "GTATGTTAGTTCTACC-1",
  "CTTGGTTAGAGTCATC-1",
  "TGCAGTTATTCCACTC-1",
  "GCCTTGCCACTTCTGA-1",
  "TCACAAGCTCGAGGCA-1",
  "GTGAAGATGAGTGCCT-1",
  "CTGCTCATTAGAGTGC-1",
  "TCTAGGCGGAAGGCAC-1",
  "GGCGGAAGTCTTGTAG-1",
  "AACTGCCTCGATAGGT-1",
  "CTCAATGACCATTAAG-1",
  "TGGCCTGACATACTCA-1",
  "GGCGTGTAATTACTTA-1",
  "AACGAATTGACCGGTT-1",
  "GGAGGTCGTAGGAGCT-1",
  "GTACAGATGACATCGC-1",
  "GTCAGCCTCCGAACCT-1",
  "GTCATCACTATTGTGC-1",
  "GTTGAGCTTCGTCCGG-1",
  "TAATTACCGGTGCGAG-1",
  "TAGGACATCAGCATGG-1",
  "TAGTGAATGTTCAGTG-1",
  "TAGTGCAGTGCGGTAG-1",
  "TAGTTCCGATAGAATC-1",
  "TCCGAATCATGTTGTC-1",
  "TGGATCATTGTAACGG-1",
  "CAGAGGCACTTAAGTC-1",
  "TGTAATCATATGAGAC-1",
  "TGTTCAGAACGGTGTA-1",
  "TGTTCTGCTCTGTCGT-1",
  "GGTCTTCCGCACATCC-1",
  "GATGAGTGAACGGCTG-1",
  "CTTAGTGGTACACGAC-1",
  "TCTGAGTCAGATCTTG-1",
  "ACTCTACTTCCGCATT-1",
  "AGACTCATGACCTGGC-1",
  "TAGACCGCGGAAGGAC-1",
  "AGCACGCTCGATCTAA-1",
  "CGCCATGAGGCGATGA-1",
  "TGAACGGATCTGGTAC-1",
  "CGAGTCATCCGCACCA-1",
  "CGGTCCATTACAACTA-1",
  "TGTCTTGGAATTGGCA-1",
  "CGTTGGTCGTCACTGT-1",
  "CGGAGCGGTAAGACGC-1",
  "CGTACATCACTGGATG-1",
  "GCCTTCTGTCACGCTC-1",
  "GCACCTAGCCGAATCT-1",
  "CAATAACTACATTAAC-1",
  "CACCATAAGCTGGCTG-1",
  "CTCCACGCTAGCTAGA-1",
  "GCATGCTGAGACAGGT-1",
  "GAGGTAACTGCGAGAC-1",
  "ACGCTGCGATTAACGA-1",
  "ATCTGGATCCTCAGAC-1",
  "CCTACACCGACATTAC-1",
  "CGCTCGCAATTGGCGA-1")


Idents(emb_29) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_29, cells_1) <- "Testis_1"
Idents(emb_29, cells_2) <- "Testis_2"
Idents(emb_29, cells_3) <- "Mesonephoros_2"
Idents(emb_29, cells_4) <- "Mesonephoros_1"
Idents(emb_29, cells_5) <- "Testis"
Idents(emb_29, cells_6) <- "Mesonephoros"

SpatialDimPlot(emb_29, pt.size.factor = 1)

de_markers <- FindMarkers(emb_29, ident.1 = "Testis")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_29_Testis_vs_ALL.csv")

de_markers <- FindMarkers(emb_29, ident.1 = "Testis", ident.2 = "Mesonephoros")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_29_Testis_ALL_vs_Mesonephoros_ALL.csv")

de_markers <- FindMarkers(emb_29, ident.1 = "Testis_1", ident.2 = "Mesonephoros_1")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_29_Testis_1_vs_Mesonephoros_1.csv")

de_markers <- FindMarkers(emb_29, ident.1 = "Testis_2", ident.2 = "Mesonephoros_2")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/emb_29_Testis_2_vs_Mesonephoros_2.csv")



