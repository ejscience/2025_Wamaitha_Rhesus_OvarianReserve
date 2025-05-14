# Defining the cell and molecular origins of the primate ovarian reserve

Sissy E. Wamaitha, Ernesto J. Rojas, Francesco Monticolo, Fei-man Hsu, Enrique Sosa, Amanda M. Mackie, Kiana Oyama, Maggie Custer, Melinda Murphy, Diana J. Laird, Jian Shu, Jon D. Hennebold, Amander T. Clark*
  
*corresponding author

[biorXiv](https://www.biorxiv.org/content/10.1101/2025.01.21.634052v1.full) DOI: 10.1101/2025.01.21.634052


## Abstract
The primate ovarian reserve is established during late fetal development and consists of quiescent primordial follicles in the ovarian cortex, each composed of granulosa cells surrounding an oocyte in dictate. As late stages of fetal development are not routinely accessible for study with human tissue, we exploited the evolutionary proximity of the rhesus macaque to investigate primate follicle formation. Similar to human prenatal ovaries, the rhesus also develops multiple types of pre-granulosa (PG) cells, with the majority of primordial follicles derived from PG2 with small variable contributions from PG1. We observed that activated medullary follicles recruit fetal theca cells to establish a two-cell system for sex-steroid hormone production prior to birth, providing a cell-based explanation for mini puberty.


## Data
Spatial transcriptomic data used in this repository and in the manuscript has been deposited in Zenodo at the following link: [insert link]. scRNA-seq data used in this manuscript has been deposited at GEO under accession number [GSE263989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi). 


## Overview
This repository details the code required to replicate the transcriptomics analyses described in the manuscript - please also see the associated Methods sections in the manuscript.


## Files
```
2025_Wamaitha_Rhesus_OvarianReserve/
  ├── README.md
├── CosMX/
│   ├── Analysis_emb020.R
│   └── Analysis_emb021.R
│   └── Session_Info
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
│   └── FilteredProbesFullRun.Rmd
│   └── Ovary_emb020.R
│   └── Ovary_emb021.R
│   └── Ovary_emb031.R
│   └── Ovary_emb032.R
│   └── Session_Info
│   └── Testis_emb019.R
│   └── Testis_emb022.R
│   └── Testis_emb029.R
├── Visium_SpaceRanger/
│   └── RhesusMacaqueProbeLists/
│       ├── filtered.probes.csv
│       ├── FRP_Human_probes_on_Macaca_mulatta.mixed_multiprobe_genes.csv
│       ├── FRP_Human_probes_on_Macaca_mulatta.multiple_hits.csv
│       ├── FRP_Human_probes_on_Macaca_mulatta.no_hit.csv
│       └── FRP_Human_probes_on_Macaca_mulatta.single_hit.csv
│       └── Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
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

To begin analysis, follow `MakingFilteredProbes.sh` to combine the single hit probe list with the subset of mixed multi-probes that are on target single hits; this generates the `filtered.probes.csv` output file.

To run spaceranger with the filtered probe list, follow `Visium_submission.sh` for the W5 and W6 samples analysed with Visium CytAssist. Follow `VisiumHD_submission.sh` for the W19 sample analysed with the Visium HD CytAssist platform.


## Visium_DEGanalysis

For the W5 and W6 samples, `FilteredProbesFullRun.Rmd` outlines the Seurat analysis workflow and sample visualisation steps used to generate Fig. 2A, 2B, S3C and S3D, and lists the required libraries and dependencies. Please note that the background H&E tissue images were converted from RGB to black and white using ImageJ.

For DEG analysis to compare gonadal regions to mesonephros or adrenal regions in the same tissue section, follow `Ovary_emb020.R` (for sample ONPRC020) or the equivalent .R file for the sample of interest (see Gonad sample key in Supplementary Materials for alternate IDs). For DEG analysis comparing gonadal regions _between_ samples, see relevant .R files in the 'Comparison' folder.


## CosMx

Analysis pipelines for the two W6 ovary samples processed on the Nanostring CosMx platform using the 1k Human Gene Panel - see manuscript Methods for additional details. Associated with Fig. 2C and S3E.


## Contact

For any questions or further information, feel free to reach out to the project lead.

---

*Updated on: May 13 2025*
