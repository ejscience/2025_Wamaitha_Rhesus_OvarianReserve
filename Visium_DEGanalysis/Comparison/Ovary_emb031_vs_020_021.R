library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

#Load EMB31
#Load_EMB31
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

#Load EMB21
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

Idents(emb_21) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_21, cells_6) <- "Ovaries"

SpatialDimPlot(emb_21, pt.size.factor = 1)

#LoadEMB_20
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

Idents(emb_20) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_20, cells_1) <- "Ovary"
SpatialDimPlot(emb_20, pt.size.factor = 1, crop = FALSE)

emb_20@meta.data[["orig.ident"]] <- "emb_20"
emb_21@meta.data[["orig.ident"]] <- "emb_21"
emb_31@meta.data[["orig.ident"]] <- "emb_31"

#Integration
embryo_list <- list(Object_name1 = emb_31, 
                    Object_name2 = emb_21,
                    Object_name3 = emb_20)

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

Idents(embryo.integrated, cells_1) <- "Ovary_31"

cells_2 <- c("CACTAGTCGACATCTA-1_2",
             "GCTTCTGTTACGAGAA-1_2",
             "GAACTCATAGTGAGTG-1_2",
             "GTCTGATAGATTATAA-1_2",
             "TCAAGGTCTGAATATG-1_2",
             "ACATGTTCCAGGACGT-1_2",
             "CCAAGCTGGCCATTGC-1_2",
             "AAGAATTCTGGCTGCA-1_2",
             "ACGGCTATGGCATACC-1_2",
             "ATGCTACTGAATAGTA-1_2",
             "TGCATTCGCTACAATG-1_2",
             "GATACTTCACTCAGCA-1_2",
             "TGCGGACTTCTGTCGC-1_2",
             "ATGTAGCCGACCAAGG-1_2",
             "GCTGAGGACAAGAAGT-1_2",
             "CTGTTCGTATGCACCA-1_2",
             "CGGATTGCGTCTACCT-1_2",
             "GTCAGTAATCCAATTA-1_2",
             "TGACCTCCTTCAGATT-1_2","AAGGCGAACTAGTAGA-1_2",
             "ACGATGCTCCAATACA-1_2",
             "AGCGGTACGCCGAGAA-1_2",
             "ATAGCTTCACAACCGT-1_2",
             "ATGACCAGACGGCCGG-1_2",
             "CAATGCGGATGCTGGT-1_2",
             "CCAGACATCGGCACGC-1_2",
             "CCGCCATACTTCACCA-1_2",
             "CCTTCCAACCTTCAGA-1_2",
             "CGACGCCAAGGTGCAA-1_2",
             "CGCCTCTAATTCGCGG-1_2",
             "CGGTAGAACGTGACTT-1_2",
             "CGGTAGATTCCTATGC-1_2",
             "CGTGGCTCCTTGTCTT-1_2",
             "CTCAGGCACGTCAGAC-1_2",
             "CTCGGAGTCTTAAGGA-1_2",
             "CTGAAGGATCTGTACG-1_2",
             "CTTACGCCATGCACAG-1_2",
             "GACGAGCTGCACTCGT-1_2",
             "GATAATCAGCTATTGT-1_2",
             "GCACTGTTGCGGCCTT-1_2",
             "GCGTACGAATCTAATC-1_2",
             "GTAGAAGTAAGCTGCG-1_2",
             "GTAGCGCTCATCAGGT-1_2",
             "GTTATGATACGATTAT-1_2",
             "GTTCCTTGAGGAACGC-1_2",
             "TAGTGACTCGCTCCGC-1_2",
             "TATGACAGCGAGAAGC-1_2",
             "TCTACAGCACGTCCGA-1_2",
             "TCTACATGAACTATAC-1_2",
             "TCTTACGGCCGCAACA-1_2",
             "TGCACGGCTGAGACTT-1_2",
             "TGGACTGAGTACCTAG-1_2",
             "TGGAGTGATGGCTCTT-1_2",
             "TGTTGGAACCTTCCGC-1_2",
             "AAGTTAAGCACTATAA-1_2"
)

Idents(embryo.integrated, cells_2) <- "Ovary_20"

cells_3 <- c("AAGCAGCGGACACGTG-1_3",
             "AATAGTATAGTTGTCG-1_3",
             "AATTGCATTGAAGAAT-1_3",
             "ACGATGTTACCTCATT-1_3",
             "ACGTTAATACAATGAC-1_3",
             "ACTCAGTGAGCCTCCA-1_3",
             "AGAAGCAGCCTAAGCT-1_3",
             "AGACCGACCGCAGACA-1_3",
             "AGACTCATTGCGGACA-1_3",
             "AGCCGGCCTATTCCAG-1_3",
             "AGCGGTGAGTGACATG-1_3",
             "AGTCGCCACTCGGCGA-1_3",
             "AGTCGGTCTAATCTGT-1_3",
             "ATACATCAAGGATTGC-1_3",
             "ATAGTTCGCAGCCGTG-1_3",
             "ATCGCAATGCACTAAC-1_3",
             "ATCTGAGGAACGGTTG-1_3",
             "ATCTGTACGCGTCTCC-1_3",
             "ATGATATCTAGTCTTA-1_3",
             "ATTACATTGAAGTATG-1_3",
             "CACGTCCACGACCTAT-1_3",
             "CAGATTGCCTGTAAGT-1_3",
             "CATCAGATACTAATCG-1_3",
             "CCTAGCCTTCATCGGC-1_3",
             "CGACAATTATTCTGTC-1_3",
             "CGAGGAAGACTTAGAT-1_3",
             "CGAGTGCGGCTCGTGG-1_3",
             "CGCTCGGTACTGCGAA-1_3",
             "CGCTTATGGTTAGGAT-1_3",
             "CGGACTGTGAGTATAG-1_3",
             "CGTCAAGGCGGAATGT-1_3",
             "CGTCCATCGTCCAATA-1_3",
             "CTACATCTACTGGCTT-1_3",
             "CTCAATACGTTAACCG-1_3",
             "CTCACGAGTAGATGAG-1_3",
             "CTCCTACATTGTATCT-1_3",
             "CTGTTCCAGAAGCTGA-1_3",
             "CTTACGGTAAGCCGAA-1_3",
             "CTTAGCTACGCAAGGT-1_3",
             "CTTCGATACTCCGACC-1_3",
             "GAATAAGGAGAACTAA-1_3",
             "GAATGACACGACTAGG-1_3",
             "GACAGTCACTTCGCTC-1_3",
             "GAGGTGCATTCGGATC-1_3",
             "GATGCCTAGAGTCTGT-1_3",
             "GCAGGCCATTCCTGCC-1_3",
             "GCCAGTAGACGATACC-1_3",
             "GCCGGTATGCTCATGC-1_3",
             "GGAATCCTTAGACTGT-1_3",
             "GGACCGCTTCCTGGCA-1_3",
             "GGTAAGGCGACGATGA-1_3",
             "GTCTAAGATCGGACTA-1_3",
             "GTCTTGCGGCGATTCA-1_3",
             "GTGCAACGACGGTGTT-1_3",
             "GTGGTTGCCGCTTCAA-1_3",
             "GTTGACTAGAGCACGT-1_3",
             "TAATAGCGCTGTACTC-1_3",
             "TAGAGGCGATTGATTC-1_3",
             "TATCCACTTGCGTATC-1_3",
             "TCAAGGTATCAACTTA-1_3",
             "TCAGTCGCCAAGGTTG-1_3",
             "TCAGTTCTCCGGACCG-1_3",
             "TCCGTGTCTTACATAA-1_3",
             "TCGATCCAGTCTGTAG-1_3",
             "TCGTCCAGCATGGCTT-1_3",
             "TCTACTGTCTCAATGA-1_3",
             "TCTTAGTCGTCTCGGC-1_3",
             "TGATAGGAGCATCAAC-1_3",
             "TGATGGCCGGCATCGC-1_3",
             "TGCAGACAATACTGTT-1_3",
             "TGGCCTGGCACTGTGG-1_3",
             "TGGTCGCCGGCTTACT-1_3",
             "TGGTGCAGACCGAAGT-1_3")
Idents(embryo.integrated, cells_3) <- "Ovary_21"

SpatialDimPlot(embryo.integrated)

embryo.integrated <- PrepSCTFindMarkers(embryo.integrated)
de_markers <- FindMarkers(embryo.integrated, ident.1 = "Ovary_31", ident.2 = "Ovary_20")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Ovary_31_vs_Ovary_20.csv")

de_markers <- FindMarkers(embryo.integrated, ident.1 = "Ovary_31", ident.2 = "Ovary_21")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Ovary_31_vs_Ovary_21.csv")

