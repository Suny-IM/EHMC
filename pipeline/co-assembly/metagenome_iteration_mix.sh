#!/bin/bash

dir=$1
threads=$2


sample_list=${dir}/sample_list
samplename=`awk -F '\t' '{print $1}' ${sample_list}`

iteration_num=0
iteration_num2=1

if [ ! -d ${dir}/Iteration/I${iteration_num2}/QC ]; then
    mkdir -p ${dir}/Iteration/I${iteration_num2}/QC/${samplename}
fi

if [ $iteration_num2 == 1 ];then
    read1=`awk -F '\t' '{print $2}' ${sample_list}`
    read2=`awk -F '\t' '{print $3}' ${sample_list}`
    ln -s ${read1} ${dir}/Iteration/I${iteration_num2}/QC/${samplename}/${samplename}_final_1.fastq
    ln -s ${read2} ${dir}/Iteration/I${iteration_num2}/QC/${samplename}/${samplename}_final_2.fastq
fi

echo "start"
date

sh /pipeline/co-assembly/assembly_binning.sh \
-i ${sample_list} \
-o ${dir}/Iteration/I${iteration_num2} \
-t ${threads} \
-s Megahit,MetaWrap,GetBbins_Binrefinement,GetBbinsUnmapped_Binrefinement_bwa \
-r yes \
-m PE

echo "end"
date

ls -lL ${dir}/Iteration/I${iteration_num2}/Binning/metawrap/*megahit/bin_refinement/B_bins_mapping/*_unmapped.1.fq | awk '{print $5,$9}' > ${dir}/Iteration/I${iteration_num2}/I${iteration_num2}_record.txt

touch ${dir}/Iteration/I${iteration_num2}/I${iteration_num2}_detail.txt

Bbin_num=`grep B_bins_unmapped.1 ${dir}/Iteration/I${iteration_num2}/I${iteration_num2}_record.txt | wc -l`


while [ ${Bbin_num} -gt 0 ] && [ ${iteration_num2} -le 2 ]
do
    iteration_num=${iteration_num2}
    iteration_num2=`expr ${iteration_num2} + 1`
    itertionDir=${dir}/Iteration/I${iteration_num2}
    mkdir -p ${itertionDir}/QC

    cd ${itertionDir}

    echo -ne "I${iteration_num2}\t${dir}/QC/I${iteration_num2}/I${iteration_num2}_final_1.fastq\t${dir}/QC/I${iteration_num2}/I${iteration_num2}_final_2.fastq\tMegahit\n" > ${dir}/Iteration/I${iteration_num2}/sample_list_I${iteration_num2}
    mkdir ${dir}/Iteration/I${iteration_num2}/QC/I${iteration_num2}
    ln -s ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_megahit/bin_refinement/B_bins_mapping/I${iteration_num}_*_unmapped.1.fq ${dir}/Iteration/I${iteration_num2}/QC/I${iteration_num2}/I${iteration_num2}_final_1.fastq
    ln -s ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_megahit/bin_refinement/B_bins_mapping/I${iteration_num}_*_unmapped.2.fq ${dir}/Iteration/I${iteration_num2}/QC/I${iteration_num2}/I${iteration_num2}_final_2.fastq



    num=1
    while read c
    do
        echo ${c} > ${itertionDir}/sample_list_I${iteration_num2}_${num}
        echo "${itertionDir}/sample_list_I${iteration_num2}_${num}" >> sampleLL
        num=`expr ${num} + 1` 
        echo ${num}
    done < ${itertionDir}/sample_list_I${iteration_num2}

    sampleCount=`cat ${itertionDir}/sampleLL | wc -l`
    useThread=`expr $threads / $sampleCount`

    while read c
    do
    {
        sh /pipeline/co-assembly/assembly_binning.sh \
        -i ${c} \
        -o ${itertionDir} \
        -t ${useThread} \
        -s Megahit,MetaWrap,GetBbins_Binrefinement,GetBbinsUnmapped_Binrefinement_bwa \
        -r yes \
        -m PE
    }&
    done < ${itertionDir}/sampleLL
    wait

    ls -lL ${itertionDir}/Binning/metawrap/*megahit/bin_refinement/B_bins_mapping/*_unmapped.1.fq | awk '{print $5,$9}' > ${itertionDir}/I${iteration_num2}_record.txt
    Bbin_num=`grep B_bins_unmapped.1 ${itertionDir}/I${iteration_num2}_record.txt | wc -l`
    
done
