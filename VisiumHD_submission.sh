
OUTPUT = /tank/data/ernesto/data/sptr/spaceranger_outs/VisiumHD

#all these are rhesus
spaceranger count \
    --id=emb017 \
    --transcriptome=/tank/data/ernesto/genomes/refdata-gex-GRCh38-2020-A \
    --probe-set=/tank/data/ernesto/data/sptr/RhMa_singlehit_mixed_single.probes.csv \
    --fastqs=/tank/data/ernesto/data/sptr/fastq/RhMa_20240530/emb017 \
    --sample=ONPRC017-D130 \
    --cytaimage=/tank/data/ernesto/data/sptr/images/HD_4samples/ONPRC017-D130_cytaimage.tiff \
    --image=/tank/data/ernesto/data/sptr/images/HD_4samples/ONPRC017-D130_image.tif \
    --slide=H1-76T7TTC \
    --area=D1 \
    --loupe-alignment=/tank/data/ernesto/data/sptr/images/HD_4samples/ONPRC017-D130_cytassist_alignment.json \
    --localcores=8 \
    --localmem=64 \
    --create-bam=true

spaceranger count \
    --id=26861 \
    --transcriptome=/tank/data/ernesto/genomes/refdata-gex-GRCh38-2020-A \
    --probe-set=/tank/data/ernesto/data/sptr/RhMa_singlehit_mixed_single.probes.csv \
    --fastqs=/tank/data/ernesto/data/sptr/fastq/RhMa_20240530/26861 \
    --sample=26861-6MPN \
    --cytaimage=/tank/data/ernesto/data/sptr/images/HD_4samples/26861-6MPN_cytaimage.tiff \
    --image=/tank/data/ernesto/data/sptr/images/HD_4samples/26861-6MPN_image.tif \
    --slide=H1-DZFRQNG \
    --area=A1 \
    --loupe-alignment=/tank/data/ernesto/data/sptr/images/HD_4samples/26861-6MPN_cytassist_alignment.json \
    --localcores=8 \
    --localmem=64 \
    --create-bam=true

######-------------
