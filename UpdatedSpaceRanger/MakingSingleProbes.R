setwd("/Users/ernestorojas/Library/Mobile Documents/com~apple~CloudDocs/Documents/UCSF/UCSF PhD/LairdMackenzie/Projects/RhMa/UpdatedSpaceRanger")

df<-read.csv(file = "RhesusMacaqueProbeLists/FRP_Human_probes_on_Macaca_mulatta.single_hit.csv")

head(df)

probe<-read.csv(file = "Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv",comment.char = "#")

head(probe)

probe$included<-FALSE


filt<-probe$probe_id%in%df$probe_id

probe[filt,4]<-TRUE

header<-read.csv(file = "Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv",nrows = 5,header = F)

write.table(x = header,file = "single.probe.csv",row.names=FALSE,col.names=FALSE,sep=",", quote = FALSE)

write.table(x = probe,file = "single.probe.csv",row.names=FALSE, sep=",", quote = FALSE,append = TRUE)
