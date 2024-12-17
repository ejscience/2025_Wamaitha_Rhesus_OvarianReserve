#first we filter for only mixed probes that are single hits and then extract the probe id. I assume you don't want probes that are predicted to hit multiple targets
grep -E -v "ENSMM.*ENSMM" RhesusMacaqueProbeLists/FRP_Human_probes_on_Macaca_mulatta.mixed_multiprobe_genes.csv |grep "ENSMM"|cut -f3 -d"|"|cut -f1 -d"," >mixedprobes_single.csv

#next we will filter for the list of single hit probes that we want to save
grep -E "ENSG" RhesusMacaqueProbeLists/FRP_Human_probes_on_Macaca_mulatta.single_hit.csv|cut -f3 -d"|"|cut -f1 -d"," > probes.txt

#We will merge the two probe sets
cat mixedprobes_single.csv probes.txt > merged_mixed_single.txt

#Now we make a new probe list with the merged list
head -n6 Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv > header.txt

#This command filters the probe set for only the ones with the identified probe ids.
grep -f merged_mixed_single.txt Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv > filtered.probes.csv

#now we make the new probeset (the header info is needed for the pipeline)
cat header.txt filtered.probes.csv > singlehit_mixed_single.probes.csv