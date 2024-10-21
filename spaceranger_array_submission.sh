#!/bin/bash
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH -o /c4/home/erojas/logs/job.%A_%a.out
#SBATCH --error=/c4/home/erojas/logs/job.%A_%a.err
#SBATCH --array=1-6

cd /c4/home/erojas

export PATH=/c4/home/erojas/software/spaceranger-2.1.0/bin:$PATH

# Array of IDs
declare -a IDs=("emb021" "emb020" "emb019" "emb022" "emb029" "emb177")

# Array of fastq paths
declare -a Fastqs=(
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/fastq/RhMa_20230726/Rh-SpTr/emb021"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/fastq/RhMa_20230726/Rh-SpTr/emb020"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/fastq/RhMa_20231103/Rh-SpTr/Rh-SpTr-d41_019"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/fastq/RhMa_20231103/Rh-SpTr/Rh-SpTr-d41_022"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/fastq/RhMa_20231103/Rh-SpTr/Rh-SpTr-d41_029"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/fastq/RhMa_20231103/Rh-SpTr/Rh-SpTr-d41_177"
)

# Array of sample names
declare -a Samples=("Rh-SpTr-d41_021" "Rh-SpTr-d41_020" "ONPRC019" "ONPRC022" "ONPRC029" "ONPRC177")

# Array of cytaimage paths
declare -a Cytaimages=(
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb020-021_images/images_emb020-021/emb021/embryo021_cytaimage.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb020-021_images/images_emb020-021/emb020/embryo020_cytaimage.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb019/emb019_cytaimage.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb022/emb022_cytaimage.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb029/emb029_cytaimage.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb177/emb177_cytaimage.tif"
)

# Array of image paths
declare -a Images=(
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb020-021_images/images_emb020-021/emb021/embryo021_slice2_image.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb020-021_images/images_emb020-021/emb020/embryo020_slice3_image.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb019/emb019_slice3_image.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb022/emb022_slice2_image.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb029/emb029_slice2_image.tif"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb177/emb177_slice3_image.tif"
)

# Array of loupe-alignment paths
declare -a LoupeAlignments=(
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb020-021_images/images_emb020-021/emb021/embryo021_cytassist_alignment.json"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb020-021_images/images_emb020-021/emb020/embryo020_cytassist_alignment.json"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb019/emb019_cytassist_alignment.json"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb022/emb022_cytassist_alignment.json"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb029/emb029_cytassist_alignment.json"
    "/c4/home/erojas/data/sptr/spatialtranscriptomics/images/RhMa_emb019-022-029-177_images/emb177/emb177_cytassist_alignment.json"
)

# Array of slide values
declare -a Slides=("V43J23-062" "V43J23-062" "V43J31-387" "V43J31-387" "V43J31-300" "V43J31-300")

# Array of Slide Areas
declare -a Areas=("A1" "D1" "D1" "A1" "D1" "A1")

# Select parameters based on SLURM_ARRAY_TASK_ID
ID=${IDs[$SLURM_ARRAY_TASK_ID-1]}
FASTQ=${Fastqs[$SLURM_ARRAY_TASK_ID-1]}
SAMPLE=${Samples[$SLURM_ARRAY_TASK_ID-1]}
CYTAIMAGE=${Cytaimages[$SLURM_ARRAY_TASK_ID-1]}
IMAGE=${Images[$SLURM_ARRAY_TASK_ID-1]}
LOUPEALIGNMENT=${LoupeAlignments[$SLURM_ARRAY_TASK_ID-1]}
SLIDE=${Slides[$SLURM_ARRAY_TASK_ID-1]}
AREA=${Areas[$SLURM_ARRAY_TASK_ID-1]}

# Run spaceranger count with selected parameters
spaceranger count \
    --id=$ID \
    --transcriptome=/c4/home/erojas/genomes/refdata-gex-GRCh38-2020-A \
    --probe-set=/c4/home/erojas/data/sptr/spatialtranscriptomics/RhMa_singlehit_mixed_single.probes.csv \
    --fastqs=$FASTQ \
    --sample=$SAMPLE \
    --cytaimage=$CYTAIMAGE \
    --image=$IMAGE \
    --slide=$SLIDE \
    --area=$AREA \
    --loupe-alignment=$LOUPEALIGNMENT \
    --localcores=8 \
    --localmem=96
