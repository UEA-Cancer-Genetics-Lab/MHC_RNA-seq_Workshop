#!/bin/bash

mkdir -p reference_files && \
    wget -nv https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz && \
    gunzip gencode.v19.pc_transcripts.fa.gz && \
    mv gencode.v19.pc_transcripts.fa reference_files/gencode.v19.pc_transcripts.fa

wget -nv https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz && \
    gunzip gencode.v19.annotation.gtf.gz && \
    mv gencode.v19.annotation.gtf reference_files/gencode.v19.annotation.gtf && \
    echo "Files in /reference_files:" && \
    ls -lh reference_files

salmon index -t reference_files/gencode.v19.pc_transcripts.fa -i reference_files/salmon_index --gencode -k 7

rm -rf references

mkdir results/fastqc
mkdir results/fastp
