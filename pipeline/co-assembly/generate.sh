file=$1
for line in `cat ${file}`;do
    array=(${line//,/ })
    sample=${array[0]}
    maindir=${array[1]}
    mkdir -p ${maindir}/QC/${sample}/remove_A_bowtie2
    rm -rf ${maindir}/QC/${sample}/remove_A_bowtie2/*
    script=${maindir}/QC/${sample}/remove_A_bowtie2/${sample}_removeA.sh
    echo "echo ${sample}_removeA" >> ${script}
    echo "date " >> ${script}
    echo "bowtie2 -p 20 \\" >> ${script}
    echo "-x ${drep_dir}/A_bins_combine/all_A_contigs_combine.fa \\" >> ${script}
    echo "-1 ${maindir}/QC/${sample}/${sample}_final_1.fastq \\" >> ${script}
    echo "-2 ${maindir}/QC/${sample}/${sample}_final_2.fastq \\" >> ${script}
    echo "-S ${maindir}/QC/${sample}/remove_A_bowtie2/${sample}_bowtie2.sam \\" >> ${script}
    echo "--un-conc ${maindir}/QC/${sample}/remove_A_bowtie2/${sample}_un_conc.fq \\" >> ${script}
    echo "--un ${maindir}/QC/${sample}/remove_A_bowtie2/${sample}_un.fq" >> ${script}
    echo "date " >> ${script}
done
date
