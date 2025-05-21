# Defining the cell and molecular origins of the primate ovarian reserve

Sissy E. Wamaitha, Ernesto J. Rojas, Francesco Monticolo, Fei-man Hsu, Enrique Sosa, Amanda M. Mackie, Kiana Oyama, Maggie Custer, Melinda Murphy, Diana J. Laird, Jian Shu, Jon D. Hennebold, Amander T. Clark*
  
*corresponding author

[biorXiv](https://www.biorxiv.org/content/10.1101/2025.01.21.634052v1.full) DOI: 10.1101/2025.01.21.634052


## Abstract
The primate ovarian reserve is established during late fetal development and consists of quiescent primordial follicles in the ovarian cortex, each composed of granulosa cells surrounding an oocyte in dictate. As late stages of fetal development are not routinely accessible for study with human tissue, we exploited the evolutionary proximity of the rhesus macaque to investigate primate follicle formation. Similar to human prenatal ovaries, the rhesus also develops multiple types of pre-granulosa (PG) cells, with the majority of primordial follicles derived from PG2 with small variable contributions from PG1. We observed that activated medullary follicles recruit fetal theca cells to establish a two-cell system for sex-steroid hormone production prior to birth, providing a cell-based explanation for mini puberty.


## Data
Rhesus macaque scRNA-seq data generated for this project is available at GEO under accession number [GSE263989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263989). Spatial transcriptomic data will be available at Zenodo [10.5281/zenodo.15477421].

The cynomolgus macaque scRNA-seq datasets included were previously published and are available at GEO ([GSE160043](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160043), [GSE194264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194264), [GSE149629](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149629)) or Zenodo ([10.5281/zenodo.6918355](https://zenodo.org/records/6918355)).


## Overview
This repository details the code required to replicate the transcriptomics analyses described in the manuscript - please also see the associated Methods sections.


## Files
```
2025_Wamaitha_Rhesus_OvarianReserve/
  ├── README.md
├── CosMX/
│   ├── Analysis_emb020.R
│   ├── Analysis_emb021.R
│   ├── Session_Info
├── Chromium/
│   ├── Analysis_Cyno.R
│   ├── Analysis_Rhesus.R
│   ├── Processing_Rhesus.R
│   ├── Session_Info
├── VisiumHD/
│   ├── Analysis_VisiumHD.R
│   ├── Session_Info
├── Visium_DEGanalysis/
│   └── Comparison/
│       ├── Ovary_emb031_vs_020_021.R
│       ├── Ovary_emb031_vs_Testis_emb029.R
│       ├── Ovary_emb032_vs_020_021.R
│       ├── Ovary_emb032_vs_031.R
│       ├── Testis_emb019_vs_022_029.R
│       ├── Testis_emb019_vs_Ovary_emb020_021.R
│       ├── Testis_emb022_vs_Ovary_emb032_031.R
│       ├── Testis_emb029_vs_Ovary_emb020_021.R
│   ├── FilteredProbesFullRun.Rmd
│   ├── Ovary_emb020.R
│   ├── Ovary_emb021.R
│   ├── Ovary_emb031.R
│   ├── Ovary_emb032.R
│   ├── Session_Info
│   ├── Testis_emb019.R
│   ├── Testis_emb022.R
│   ├── Testis_emb029.R
├── Visium_SpaceRanger/
│   └── RhesusMacaqueProbeLists/
│       ├── filtered.probes.csv
│       ├── FRP_Human_probes_on_Macaca_mulatta.mixed_multiprobe_genes.csv
│       ├── FRP_Human_probes_on_Macaca_mulatta.multiple_hits.csv
│       ├── FRP_Human_probes_on_Macaca_mulatta.no_hit.csv
│       ├── FRP_Human_probes_on_Macaca_mulatta.single_hit.csv
│       ├── Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
│   ├── MakingFilteredProbes.sh
│   ├── Session_Info
│   ├── VisiumHD_submission.sh
│   ├── Visium_submission.sh
```


## Methods

## Visium_SpaceRanger

The `RhesusMacaqueProbeLists` folder contains the results of the in silico analyisis provided by 10x Genomics, which compares the 10x Visium Human Panel probe set to predicted targets in the rhesus macaque:
  
- Whole probe set: `Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv`
- Single hit (on target): `FRP_Human_probes_on_Macaca_mulatta.single_hit.csv`
- Mixed multi-probe hits (some on target single hits, some off target): `FRP_Human_probes_on_Macaca_mulatta.mixed_multiprobe_genes.csv`
- Multiple hits (none on target): `FRP_Human_probes_on_Macaca_mulatta.multiple_hits.csv`
- No hits: `FRP_Human_probes_on_Macaca_mulatta.no_hit.csv`

To begin analysis, we followed `MakingFilteredProbes.sh` to combine the single hit probe list with the subset of mixed multi-probes that are on target single hits; this generates the `filtered.probes.csv` output file.

To run spaceranger with the filtered probe list, we followed `Visium_submission.sh` for the W5 and W6 samples analysed with Visium CytAssist, and `VisiumHD_submission.sh` for the W19 sample analysed with the Visium HD CytAssist platform.


## Visium_DEGanalysis

For the W5 and W6 samples, `FilteredProbesFullRun.Rmd` outlines the Seurat analysis workflow and sample visualisation steps used to generate Fig. 2A, 2B, S3C and S3D, and lists the required libraries and dependencies. Please note that the background H&E tissue images were converted from RGB to black and white using ImageJ.

DEG analysis to compare gonadal regions to mesonephros or adrenal regions in the same tissue section is outlined in `Ovary_emb020.R` (for sample ONPRC020) or the equivalent .R file for the sample of interest (see Gonad sample key in Supplementary Materials for alternate IDs). For DEG analysis comparing gonadal regions _between_ samples, see relevant .R files in the 'Comparison' folder.


## CosMx

Analysis pipelines for the two W6 ovary samples processed on the Nanostring CosMx platform using the 1k Human Gene Panel - see manuscript Methods for additional details.  Associated with Fig. 2C and S3E.


## Chromium

Pre-processing steps for the W8 - W19 rhesus macaque scRNA seq samples are detailed in `Processing_Rhesus.R`; we followed the steps in `Analysis_Rhesus.R` for the Seurat workflow and downstream analysis. 

For the analysis of the published cynomolgus macaque datasets, we followed the steps in `Analysis_Cyno.R` to combine, process and integrate the datasets (Harmony integration) as well as downstream analysis. Associated with Fig. S4F-G and S5C-E.


## VisiumHD

Analysis pipeline for the W19 ovary sample processed on the Visium HD CytAssist platform - see manuscript Methods for additional details.  Associated with Fig. 4C, 5A, S6E and S6F.


## Contact

For any questions or further information, reach out to the project lead.

---
  
  *Updated on: May 21 2025*
  
