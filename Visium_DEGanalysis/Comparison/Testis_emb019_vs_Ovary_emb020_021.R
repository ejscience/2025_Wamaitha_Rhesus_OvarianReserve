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

####OVARY 20#####
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

Idents(emb_20) <- "ident"  # Replace "ident" with the actual identity class in your Seurat object
Idents(emb_20, cells_1) <- "Ovary"
Idents(emb_20, cells_2) <- "Mesonephoros"
SpatialDimPlot(emb_20, pt.size.factor = 1, crop = FALSE)

####OVARY 21#####
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
Idents(emb_21, cells_6) <- "Ovaries"
Idents(emb_21, cells_7) <- "Mesonephoros"

SpatialDimPlot(emb_21, pt.size.factor = 1)


emb_19@meta.data[["orig.ident"]] <- "emb_19"
emb_20@meta.data[["orig.ident"]] <- "emb_20"
emb_21@meta.data[["orig.ident"]] <- "emb_21"

#Integration
embryo_list <- list(Object_name1 = emb_19, 
                    Object_name2 = emb_20,
                    Object_name2 = emb_21)

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


cells_2 <- c("AAGCAGCGGACACGTG-1_2",
             "AATAGTATAGTTGTCG-1_2",
             "AATTGCATTGAAGAAT-1_2",
             "ACGATGTTACCTCATT-1_2",
             "ACGTTAATACAATGAC-1_2",
             "ACTCAGTGAGCCTCCA-1_2",
             "AGAAGCAGCCTAAGCT-1_2",
             "AGACCGACCGCAGACA-1_2",
             "AGACTCATTGCGGACA-1_2",
             "AGCCGGCCTATTCCAG-1_2",
             "AGCGGTGAGTGACATG-1_2",
             "AGTCGCCACTCGGCGA-1_2",
             "AGTCGGTCTAATCTGT-1_2",
             "ATACATCAAGGATTGC-1_2",
             "ATAGTTCGCAGCCGTG-1_2",
             "ATCGCAATGCACTAAC-1_2",
             "ATCTGAGGAACGGTTG-1_2",
             "ATCTGTACGCGTCTCC-1_2",
             "ATGATATCTAGTCTTA-1_2",
             "ATTACATTGAAGTATG-1_2",
             "CACGTCCACGACCTAT-1_2",
             "CAGATTGCCTGTAAGT-1_2",
             "CATCAGATACTAATCG-1_2",
             "CCTAGCCTTCATCGGC-1_2",
             "CGACAATTATTCTGTC-1_2",
             "CGAGGAAGACTTAGAT-1_2",
             "CGAGTGCGGCTCGTGG-1_2",
             "CGCTCGGTACTGCGAA-1_2",
             "CGCTTATGGTTAGGAT-1_2",
             "CGGACTGTGAGTATAG-1_2",
             "CGTCAAGGCGGAATGT-1_2",
             "CGTCCATCGTCCAATA-1_2",
             "CTACATCTACTGGCTT-1_2",
             "CTCAATACGTTAACCG-1_2",
             "CTCACGAGTAGATGAG-1_2",
             "CTCCTACATTGTATCT-1_2",
             "CTGTTCCAGAAGCTGA-1_2",
             "CTTACGGTAAGCCGAA-1_2",
             "CTTAGCTACGCAAGGT-1_2",
             "CTTCGATACTCCGACC-1_2",
             "GAATAAGGAGAACTAA-1_2",
             "GAATGACACGACTAGG-1_2",
             "GACAGTCACTTCGCTC-1_2",
             "GAGGTGCATTCGGATC-1_2",
             "GATGCCTAGAGTCTGT-1_2",
             "GCAGGCCATTCCTGCC-1_2",
             "GCCAGTAGACGATACC-1_2",
             "GCCGGTATGCTCATGC-1_2",
             "GGAATCCTTAGACTGT-1_2",
             "GGACCGCTTCCTGGCA-1_2",
             "GGTAAGGCGACGATGA-1_2",
             "GTCTAAGATCGGACTA-1_2",
             "GTCTTGCGGCGATTCA-1_2",
             "GTGCAACGACGGTGTT-1_2",
             "GTGGTTGCCGCTTCAA-1_2",
             "GTTGACTAGAGCACGT-1_2",
             "TAATAGCGCTGTACTC-1_2",
             "TAGAGGCGATTGATTC-1_2",
             "TATCCACTTGCGTATC-1_2",
             "TCAAGGTATCAACTTA-1_2",
             "TCAGTCGCCAAGGTTG-1_2",
             "TCAGTTCTCCGGACCG-1_2",
             "TCCGTGTCTTACATAA-1_2",
             "TCGATCCAGTCTGTAG-1_2",
             "TCGTCCAGCATGGCTT-1_2",
             "TCTACTGTCTCAATGA-1_2",
             "TCTTAGTCGTCTCGGC-1_2",
             "TGATAGGAGCATCAAC-1_2",
             "TGATGGCCGGCATCGC-1_2",
             "TGCAGACAATACTGTT-1_2",
             "TGGCCTGGCACTGTGG-1_2",
             "TGGTCGCCGGCTTACT-1_2",
             "TGGTGCAGACCGAAGT-1_2")

Idents(embryo.integrated, cells_2) <- "Ovary_20"

cells_3 <- c("AAGAATTCTGGCTGCA-1_3",
             "AAGGCGAACTAGTAGA-1_3",
             "AAGTTAAGCACTATAA-1_3",
             "ACATGTTCCAGGACGT-1_3",
             "ACGATGCTCCAATACA-1_3",
             "ACGGCTATGGCATACC-1_3",
             "AGCGGTACGCCGAGAA-1_3",
             "ATAGCTTCACAACCGT-1_3",
             "ATGACCAGACGGCCGG-1_3",
             "ATGCTACTGAATAGTA-1_3",
             "ATGTAGCCGACCAAGG-1_3",
             "CAATGCGGATGCTGGT-1_3",
             "CACTAGTCGACATCTA-1_3",
             "CCAAGCTGGCCATTGC-1_3",
             "CCAGACATCGGCACGC-1_3",
             "CCGCCATACTTCACCA-1_3",
             "CCTTCCAACCTTCAGA-1_3",
             "CGACGCCAAGGTGCAA-1_3",
             "CGCCTCTAATTCGCGG-1_3",
             "CGGATTGCGTCTACCT-1_3",
             "CGGTAGAACGTGACTT-1_3",
             "CGGTAGATTCCTATGC-1_3",
             "CGTGGCTCCTTGTCTT-1_3",
             "CTCAGGCACGTCAGAC-1_3",
             "CTCGGAGTCTTAAGGA-1_3",
             "CTGAAGGATCTGTACG-1_3",
             "CTGTTCGTATGCACCA-1_3",
             "CTTACGCCATGCACAG-1_3",
             "GAACTCATAGTGAGTG-1_3",
             "GACGAGCTGCACTCGT-1_3",
             "GATAATCAGCTATTGT-1_3",
             "GATACTTCACTCAGCA-1_3",
             "GCACTGTTGCGGCCTT-1_3",
             "GCGTACGAATCTAATC-1_3",
             "GCTGAGGACAAGAAGT-1_3",
             "GCTTCTGTTACGAGAA-1_3",
             "GTAGAAGTAAGCTGCG-1_3",
             "GTAGCGCTCATCAGGT-1_3",
             "GTCTGATAGATTATAA-1_3",
             "GTTATGATACGATTAT-1_3",
             "GTTCCTTGAGGAACGC-1_3",
             "TAGTGACTCGCTCCGC-1_3",
             "TATGACAGCGAGAAGC-1_3",
             "TCAAGGTCTGAATATG-1_3",
             "TCTACAGCACGTCCGA-1_3",
             "TCTACATGAACTATAC-1_3",
             "TCTTACGGCCGCAACA-1_3",
             "TGACCTCCTTCAGATT-1_3",
             "TGCACGGCTGAGACTT-1_3",
             "TGCATTCGCTACAATG-1_3",
             "TGCGGACTTCTGTCGC-1_3",
             "TGGACTGAGTACCTAG-1_3",
             "TGGAGTGATGGCTCTT-1_3",
             "TGTTGGAACCTTCCGC-1_3",
             "GTCAGTAATCCAATTA-1_3"
)

Idents(embryo.integrated, cells_3) <- "Ovary_21"

SpatialDimPlot(embryo.integrated)

embryo.integrated <- PrepSCTFindMarkers(embryo.integrated)
de_markers <- FindMarkers(embryo.integrated, ident.1 = "Testis_19", ident.2 = "Ovary_20")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Testis19_vs_Ovary_20.csv")

de_markers <- FindMarkers(embryo.integrated, ident.1 = "Testis_19", ident.2 = "Ovary_21")
write.csv(de_markers, file = "~/Desktop/NEW_Visium_Rm/Results/embryo.integrated_Testis19_vs_Ovary_21.csv")
