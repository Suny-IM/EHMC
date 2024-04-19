#!/bin/bash

outputdir=$1
sample=$2
bbmap_path=$3
ref=$4


##### Step 1 : QC for pair end reads #####

read1outpaired=${outputdir}/QC/${sample}/${sample}_trim_1.fq
read2outpaired=${outputdir}/QC/${sample}/${sample}_trim_2.fq 

sh bbduk.sh \
in=$rawread1 \
in2=$rawread2 \
out=$read1outpaired \
out2=$read2outpaired  \
ref=${bbmap_path}/adapters.fa \
ktrim=r \
k=23 \
mink=11 \
hdist=1 \
tbo \
trimq=15 \
qtrim=rl \
maq=20 \
maxns=0 \
minlength=40


##### Step 1 : QC for single end reads #####

readoutpaired=${outputdir}/QC/${sample}/${sample}_trim.fq

sh bbduk.sh \
in=$rawread \
out=$readoutpaired \
ref=${bbmap_path}/resources/adapters.fa \
ktrim=r \
k=23 \
mink=11 \
hdist=1 \
tbo \
trimq=15 \
maq=20 \
maxns=0 \
minlength=40


##### Step 2: filter out contaminations #####

##Filter potential host sequence

##Filter potential host sequence for pair end reads

read1outpaired=${outputdir}/QC/${sample}/${sample}_trim_1.fq
read2outpaired=${outputdir}/QC/${sample}/${sample}_trim_2.fq
bowtie2 -p ${threads} \
-x ${Ref} \
-1 ${read1outpaired} \
-2 ${read2outpaired} \
-S $sam \
--un-conc ${outputdir}/QC/${sample}/${sample}_un_conc.fq \
--un ${outputdir}/QC/${sample}/${sample}_un.fq

mv ${outputdir}/QC/${sample}/${sample}_un_conc.1.fq ${outputdir}/QC/${sample}/${sample}_final_1.fastq
mv ${outputdir}/QC/${sample}/${sample}_un_conc.2.fq ${outputdir}/QC/${sample}/${sample}_final_2.fastq


##Filter potential host sequence for single end reads

read1outpaired=${outputdir}/QC/${sample}/${sample}_trim_1.fq
read2outpaired=${outputdir}/QC/${sample}/${sample}_trim_2.fq
bowtie2 -p ${threads} \
-x ${Ref} \
-1 ${read1outpaired} \
-2 ${read2outpaired} \
-S $sam \
--un-conc ${outputdir}/QC/${sample}/${sample}_un_conc.fq \
--un ${outputdir}/QC/${sample}/${sample}_un.fq

mv ${outputdir}/QC/${sample}/${sample}_un_conc.1.fq ${outputdir}/QC/${sample}/${sample}_final_1.fastq
mv ${outputdir}/QC/${sample}/${sample}_un_conc.2.fq ${outputdir}/QC/${sample}/${sample}_final_2.fastq

readoutpaired=${outputdir}/QC/${sample}/${sample}_trim.fq
bowtie2 -p ${threads} \
-x ${ref} \
-U ${readoutpaired} \
-S $sam \
--un ${outputdir}/QC/${sample}/${sample}_un.fq

mv ${outputdir}/QC/${sample}/${sample}_un.fq ${outputdir}/QC/${sample}/${sample}_final.fastq


##Fast QC report for paired end reads

fastqc  \
-t ${threads} \
${outputdir}/QC/${sample}/${sample}_final_1.fastq \
${outputdir}/QC/${sample}/${sample}_final_2.fastq \
-o  ${outputdir}/QC/${sample}/clean_qc

##Fast QC report for single end reads

fastqc  \
-t ${threads} \
${outputdir}/QC/${sample}/${sample}_final.fastq \
-o  ${outputdir}/QC/${sample}/clean_qc