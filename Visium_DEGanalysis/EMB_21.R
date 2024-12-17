library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

#Load
emb_21 <- Load10X_Spatial('~/Desktop/New_Visium_Rm/emb021_outs/outs/',filename = "filtered_feature_bc_matrix.h5", slice = "slice2")

#QC
VlnPlot(emb_21, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, ncol = 2) + NoLegend()
SpatialFeaturePlot(emb_21, features = c("nCount_Spatial", "nFeature_Spatial"))
emb_21 <- emb_21[, emb_21$nFeature_Spatial > 100 & emb_21$nCount_Spatial > 100]
SpatialFeaturePlot(emb_21, features = c("nCount_Spatial", "nFeature_Spatial"))
dim(emb_21)

#SCT trasnform
emb_21 <- SCTransform(emb_21, assay = "Spatial", vst.flavor = "v2", verbose = FALSE)

emb_21 <- RunPCA(emb_21, assay = "SCT", verbose = FALSE)
emb_21 <- FindNeighbors(emb_21, reduction = "pca", dims = 1:30)
emb_21 <- FindClusters(emb_21, verbose = FALSE)
emb_21 <- RunUMAP(emb_21, reduction = "pca", dims = 1:30)

SpatialDimPlot(emb_21, label = TRUE, label.size = 3, pt.size.factor = 1)

LinkedDimPlot(emb_21)

cells_2 <- c(
  "ATACGGCGACTGGTCT-1",
  "AGGCTGCAGCTCAGGC-1",
  "GTTCATACCTTGGTTA-1",
  "CGACTTCACTCGGCAT-1",
  "ATCTGGATCCTCAGAC-1",
  "TAGTGCAGTGCGGTAG-1",
  "ACGCTGCGATTAACGA-1",
  "GGAGGTCGTAGGAGCT-1",
  "CGCTCGCAATTGGCGA-1",
  "GCTGACTTGTCTCAAT-1",
  "TCACCATAGATTGTTA-1",
  "CCTGAAGCGTATAGCA-1",
  "CAATAACTACATTAAC-1",
  "TGTTCTGCTCTGTCGT-1",
  "CACCATAAGCTGGCTG-1",
  "TCGGACATTCTATGAG-1",
  "TCGATCTAATGTTAGT-1",
  "CCTGTCCTCAGTTGTA-1",
  "GAGGTAACTGCGAGAC-1",
  "TGGATCATTGTAACGG-1",
  "TGTAATCATATGAGAC-1",
  "TAGTTCCGATAGAATC-1",
  "ATGCCGTTGTGAAGCG-1",
  "GTTGAGCTTCGTCCGG-1",
  "AGCATGCACCGCGTCA-1"
)

cells_1 <- c("AACCTTAGCGTCCGGT-1",
             "GTAGCAGAGCAATAGC-1",
             "GGCCGCCAACTATTGC-1",
             "TAATTGATTCTGTCGC-1",
             "TATTAAGCAGCGCCGG-1",
             "TGTATGATTGTAGGTG-1",
             "CGTGAGCAATGGACAT-1",
             "TCTACATGAATCTGAT-1",
             "GGAAGTCTTGTCGCAT-1",
             "CCGCGTATCGGTTGCAG-1",
             "GAATTAACCTTCAATT-1",
             "AGACCGGAGCAACGCT-1",
             "TACGGATCCTTGCGAT-1",
             "ATAATTAGTTCTACGG-1",
             "AAGTACTAAGGTTACG-1",
             "ACTTAAGCCAGTAGGA-1",
             "CTATAACTCCTGGAGG-1",
             "GATCTTGTCATGATCT-1",
             "TCGAATCTTGCGTTAC-1",
             "GCGATGTACTGCTATT-1",
             "GTTGCCACGACAAGAA-1",
             "CGTCCATCCTAGCATC-1",
             "GATGACTCTGAATCTT-1",
             "GGCCGTCCAGTAGGAC-1",
             "CGGTCTAGATCGAATG-1",
             "CCGCGTATCGTTGCAG-1",
             "ATAGTTAAGCCTGCCG-1",
             "GATATGCTCACGTTAC-1",
             "TCAACATGGATACAAC-1",
             "AATTGTAGAATCCTAG-1",
             "TATGTTGCAGTCTAGG-1",
             "TGCGTGTCAGCCAGTG-1",
             "AATTCTTGGATACGGA-1",
             "AAGACAATCTCCTGCA-1",
             "GTCGCACGCAGCCTGC-1",
             "TCATAGGTAGTGACGA-1")

cells_4 <- c(
  "GTACTGCATAAGGTTA-1",
  "CGGTCTTGCTATTCTA-1",
  "CTGGCGGACAACGCAC-1",
  "TCATGAGGTTGATATG-1",
  "TGCGTATCTTGTACCG-1",
  "GCAAGAATGTTCCACC-1",
  "CTCACGATCTCCGGAG-1",
  "GGTTCTCCTATCAGTT-1",
  "ACATCGTAGACAGGTT-1",
  "CGACGCCAAGGTGCAA-1",
  "GACGGTACCTGGCCGT-1",
  "CCGGCGCTAAGGATGA-1",
  "CAACTGACAAGACAAG-1",
  "GTCTATTGCGAGTTAC-1",
  "GCACGCAACCATACCA-1",
  "CGGTCAACGGCGGATT-1",
  "ATGCATGTCACTCGTG-1",
  "CTTAGTAAGCGTGGCA-1",
  "TGATCCGCACTCCTTA-1",
  "CGTAGCCAAGTCCGTC-1",
  "AATTCTGGATACGGA-1",
  "AAGACCAATCTCCTGCA-1",
  "CTCCTGTCGAAGCCAA-1",
  "ATAGTAAGCCTGCCG-1",
  "TCAGTCCAACGGAAGT-1",
  "CAGCCTCGTGGCTCAC-1")


cells_3 <- c(
  "AAGGCGAACTAGTAGA-1",
  "ACGATGCTCCAATACA-1",
  "AGCGGTACGCCGAGAA-1",
  "ATAGCTTCACAACCGT-1",
  "ATGACCAGACGGCCGG-1",
  "CAATGCGGATGCTGGT-1",
  "CCAGACATCGGCACGC-1",
  "CCGCCATACTTCACCA-1",
  "CCTTCCAACCTTCAGA-1",
  "CGACGCCAAGGTGCAA-1",
  "CGCCTCTAATTCGCGG-1",
  "CGGTAGAACGTGACTT-1",
  "CGGTAGATTCCTATGC-1",
  "CGTGGCTCCTTGTCTT-1",
  "CTCAGGCACGTCAGAC-1",
  "CTCGGAGTCTTAAGGA-1",
  "CTGAAGGATCTGTACG-1",
  "CTTACGCCATGCACAG-1",
  "GACGAGCTGCACTCGT-1",
  "GATAATCAGCTATTGT-1",
  "GCACTGTTGCGGCCTT-1",
  "GCGTACGAATCTAATC-1",
  "GTAGAAGTAAGCTGCG-1",
  "GTAGCGCTCATCAGGT-1",
  "GTTATGATACGATTAT-1",
  "GTTCCTTGAGGAACGC-1",
  "TAGTGACTCGCTCCGC-1",
  "TATGACAGCGAGAAGC-1",
  "TCTACAGCACGTCCGA-1",
  "TCTACATGAACTATAC-1",
  "TCTTACGGCCGCAACA-1",
  "TGCACGGCTGAGACTT-1",
  "TGGACTGAGTACCTAG-1",
  "TGGAGTGATGGCTCTT-1",
  "TGTTGGAACCTTCCGC-1",
  "AAGTTAAGCACTATAA-1"
)

cells_5 <- c("CACTAGTCGACATCTA-1",
             "GCTTCTGTTACGAGAA-1",
             "GAACTCATAGTGAGTG-1",
             "GTCTGATAGATTATAA-1",
             "TCAAGGTCTGAATATG-1",
             "ACATGTTCCAGGACGT-1",
             "CCAAGCTGGCCATTGC-1",
             "AAGAATTCTGGCTGCA-1",
             "ACGGCTATGGCATACC-1",
             "ATGCTACTGAATAGTA-1",
             "TGCATTCGCTACAATG-1",
             "GATACTTCACTCAGCA-1",
             "TGCGGACTTCTGTCGC-1",
             "ATGTAGCCGACCAAGG-1",
             "GCTGAGGACAAGAAGT-1",
             "CTGTTCGTATGCACCA-1",
             "CGGATTGCGTCTACCT-1",
             "GTCAGTAATCCAATTA-1",
             "TGACCTCCTTCAGATT-1"
)


cells_6 <- c("CACTAGTCGACATCTA-1",
             "GCTTCTGTTACGAGAA-1",
             "GAACTCATAGTGAGTG-1",
             "GTCTGATAGATTATAA-1",
             "TCAAGGTCTGAATATG-1",
             "ACATGTTCCAGGACGT-1",
             "CCAAGCTGGCCATTGC-1",
             "AAGAATTCTGGCTGCA-1",
             "ACGGCTATGGCATACC-1",
             "ATGCTACTGAATAGTA-1",
             "TGCATTCGCTACAATG-1",
             "GATACTTCACTCAGCA-1",
             "TGCGGACTTCTGTCGC-1",
             "ATGTAGCCGACCAAGG-1",
             "GCTGAGGACAAGAAGT-1",
             "CTGTTCGTATGCACCA-1",
             "CGGATTGCGTCTACCT-1",
             "GTCAGTAATCCAATTA-1",
             "TGACCTCCTTCAGATT-1","AAGGCGAACTAGTAGA-1",
             "ACGATGCTCCAATACA-1",
             "AGCGGTACGCCGAGAA-1",
             "ATAGCTTCACAACCGT-1",
             "ATGACCAGACGGCCGG-1",
             "CAATGCGGATGCTGGT-1",
             "CCAGACATCGGCACGC-1",
             "CCGCCATACTTCACCA-1",
             "CCTTCCAACCTTCAGA-1",
             "CGACGCCAAGGTGCAA-1",
             "CGCCTCTAATTCGCGG-1",
             "CGGTAGAACGTGACTT-1",
             "CGGTAGATTCCTATGC-1",
             "CGTGGCTCCTTGTCTT-1",
             "CTCAGGCACGTCAGAC-1",
             "CTCGGAGTCTTAAGGA-1",
             "CTGAAGGATCTGTACG-1",
             "CTTACGCCATGCACAG-1",
             "GACGAGCTGCACTCGT-1",
             "GATAATCAGCTATTGT-1",
             "GCACTGTTGCGGCCTT-1",
             "GCGTACGAATCTAATC-1",
             "GTAGAAGTAAGCTGCG-1",
             "GTAGCGCTCATCAGGT-1",
             "GTTATGATACGATTAT-1",
             "GTTCCTTGAGGAACGC-1",
             "TAGTGACTCGCTCCGC-1",
             "TATGACAGCGAGAAGC-1",
             "TCTACAGCACGTCCGA-1",
             "TCTACATGAACTATAC-1",
             "TCTTACGGCCGCAACA-1",
             "TGCACGGCTGAGACTT-1",
             "TGGACTGAGTACCTAG-1",
             "TGGAGTGATGGCTCTT-1",
             "TGTTGGAACCTTCCGC-1",
             "AAGTTAAGCACTATAA-1"
)

cells_7 <- c(
  "GTACTGCATAAGGTTA-1",
  "CGGTCTTGCTATTCTA-1",
  "CTGGCGGACAACGCAC-1",
  "TCATGAGGTTGATATG-1",
  "TGCGTATCTTGTACCG-1",
  "GCAAGAATGTTCCACC-1",
  "CTCACGATCTCCGGAG-1",
  "GGTTCTCCTATCAGTT-1",
  "ACATCGTAGACAGGTT-1",
  "CGACGCCAAGGTGCAA-1",
  "GACGGTACCTGGCCGT-1",
  "CCGGCGCTAAGGATGA-1",
  "CAACTGACAAGACAAG-1",
  "GTCTATTGCGAGTTAC-1",
  "GCACGCAACCATACCA-1",
  "CGGTCAACGGCGGATT-1",
  "ATGCATGTCACTCGTG-1",
  "CTTAGTAAGCGTGGCA-1",
  "TGATCCGCACTCCTTA-1",
  "CGTAGCCAAGTCCGTC-1",
  "AATTCTGGATACGGA-1",
  "AAGACCAATCTCCTGCA-1",
  "CTCCTGTCGAAGCCAA-1",
  "ATAGTAAGCCTGCCG-1",
  "TCAGTCCAACGGAAGT-1",
  "CAGCCTCGTGGCTCAC-1",
  "AACCTTAGCGTCCGGT-1",
  "GTAGCAGAGCAATAGC-1",
  "GGCCGCCAACTATTGC-1",
  "TAATTGATTCTGTCGC-1",
  "TATTAAGCAGCGCCGG-1",
  "TGTATGATTGTAGGTG-1",
  "CGTGAGCAATGGACAT-1",
  "TCTACATGAATCTGAT-1",
  "GGAAGTCTTGTCGCAT-1",
  "CCGCGTATCGGTTGCAG-1",
  "GAATTAACCTTCAATT-1",
  "AGACCGGAGCAACGCT-1",
  "TACGGATCCTTGCGAT-1",
  "ATAATTAGTTCTACGG-1",
  "AAGTACTAAGGTTACG-1",
  "ACTTAAGCCAGTAGGA-1",
  "CTATAACTCCTGGAGG-1",
  "GATCTTGTCATGATCT-1",
  "TCGAATCTTGCGTTAC-1",
  "GCGATGTACTGCTATT-1",
  "GTTGCCACGACAAGAA-1",
  "CGTCCATCCTAGCATC-1",
  "GATGACTCTGAATCTT-1",
  "GGCCGTCCAGTAGGAC-1",
  "CGGTCTAGATCGAATG-1",
  "CCGCGTATCGTTGCAG-1",
  "ATAGTTAAGCCTGCCG-1",
  "GATATGCTCACGTTAC-1",
  "TCAACATGGATACAAC-1",
  "AATTGTAGAATCCTAG-1",
  "TATGTTGCAGTCTAGG-1",
  "TGCGTGTCAGCCAGTG-1",
  "AATTCTTGGATACGGA-1",
  "AAGACAATCTCCTGCA-1",
  "GTCGCACGCAGCCTGC-1",
  "TCATAGGTAGTGACGA-1")

Idents(emb_21) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_21, cells_2) <- "Adrenal_gland_2"
Idents(emb_21, cells_1) <- "Mesonephoros_1"
Idents(emb_21, cells_3) <- "Ovary_2"
Idents(emb_21, cells_4) <- "Mesonephoros_2"
Idents(emb_21, cells_5) <- "Ovary_1"
Idents(emb_21, cells_6) <- "Ovaries"
Idents(emb_21, cells_7) <- "Mesonephoros"

SpatialDimPlot(emb_21, pt.size.factor = 1)

de_markers <- FindMarkers(emb_21, ident.1 = "Ovaries")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Emb_21_Ovaries_vs_ALL.csv")

de_markers <- FindMarkers(emb_21, ident.1 = "Ovaries", ident.2 = "Mesonephoros")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Emb_21_Ovary_ALL_vs_Mesonephoros_ALL.csv")

de_markers <- FindMarkers(emb_21, ident.1 = "Ovaries", ident.2 = "Adrenal_gland")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Emb_21_Ovary_ALL_vs_Adrenal_gland_ALL.csv")

de_markers <- FindMarkers(emb_21, ident.1 = "Ovary_1", ident.2 = "Mesonephoros_1")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Emb_21_Ovary_1_vs_Mesonephoros_1.csv")

de_markers <- FindMarkers(emb_21, ident.1 = "Ovary_2", ident.2 = "Mesonephoros_2")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Emb_21_Ovary_2_vs_Mesonephoros_2.csv")

de_markers <- FindMarkers(emb_21, ident.1 = "Ovary_2", ident.2 = "Adrenal_gland_2")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/Emb_21_Ovary_2_vs_Adrenal_gland_2.csv")


