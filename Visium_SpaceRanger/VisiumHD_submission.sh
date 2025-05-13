
OUTPUT = /tank/data/ernesto/data/sptr/spaceranger_outs/VisiumHD

#rhesus
spaceranger count \
    --id=emb017 \
    --transcriptome=/tank/data/ernesto/genomes/refdata-gex-GRCh38-2020-A \
    --probe-set=/tank/data/ernesto/data/sptr/RhMa_singlehit_mixed_single.probes.csv \
    --fastqs=/tank/data/ernesto/data/sptr/fastq/RhMa_20240530/emb017 \
    --sample=ONPRC017 \
    --cytaimage=/tank/data/ernesto/data/sptr/images/HD_4samples/ONPRC017_cytaimage.tiff \
    --image=/tank/data/ernesto/data/sptr/images/HD_4samples/ONPRC017_image.tif \
    --slide=H1-76T7TTC \
    --area=D1 \
    --loupe-alignment=/tank/data/ernesto/data/sptr/images/HD_4samples/ONPRC017_cytassist_alignment.json \
    --localcores=8 \
    --localmem=64 \
    --create-bam=true

######-------------
