#!/bin/bash

while getopts ":i:o:t:s:d:p:r:m:n:" opt ;do
    case $opt in
        i) inputfile=$OPTARG;;
        o) outputdir=$OPTARG;;
        t) threads=$OPTARG;;
        s) software=$OPTARG;;
        r) bin_refinement=$OPTARG;;
        m) mode=$OPTARG;;
    esac
done

echo $mode

if [ -n "${software}" ];then
    array=(${software//,/ })
fi


if [ ! -n "${threads}" ];then
    threads=30
fi


if [ ! -d ${outputdir} ]; then
        mkdir -p ${outputdir}
fi


if [ ! -d ${outputdir}/work_dir ]; then
        mkdir -p ${outputdir}/work_dir
fi




##########load_configFile


if [ $mode == "PE" ]; then
    while read line
    do
            tmp=($line)
            name=(${name[@]} ${tmp[0]})
            pe1=(${pe1[@]} ${tmp[1]})
            pe2=(${pe2[@]} ${tmp[2]})
            startpoints=(${startpoints[@]} ${tmp[3]})
    done < $inputfile
elif [ $mode == "SE" ]; then
    while read line
    do
            tmp=($line)
            name=(${name[@]} ${tmp[0]})
            se=(${se[@]} ${tmp[1]})
            startpoints=(${startpoints[@]} ${tmp[2]})
    done < $inputfile
else
    echo "error! chose one for mode from PE or SE"
fi




function contains(){
    local n=$#
    local value=${!n}
    for ((i=1;i < $#;i++)) {
        if [ "${!i}" == "${value}" ]; then
            echo "y"
            return 0
        fi
    }
    echo "n"
    return 1
}


function Megahit(){

    if [ -d ${outputdir}/Assembly/${sample} ]; then
        rm -rf ${outputdir}/Assembly/${sample}
    fi
    mkdir -p ${outputdir}/Assembly/${sample}

    MegahitstartTime=`date +'%Y-%m-%d %H:%M:%S'`
    MegahitstartTime_s=$(date --date="$MegahitstartTime" +%s)

    if [ $mode == "PE" ]; then
        megahit \
        -1 ${outputdir}/QC/${sample}/${sample}_final_1.fastq \
        -2 ${outputdir}/QC/${sample}/${sample}_final_2.fastq \
        -t ${threads} \
        --k-list 21,29,39,59,79,99,119,141 \
        --out-prefix ${sample} \
        -o ${outputdir}/Assembly/${sample}/megahit_out
    elif [ $mode == "SE" ]; then
        megahit \
        -r ${outputdir}/QC/${sample}/${sample}_final.fastq \
        -t ${threads} \
        --k-list 21,29,39,59,79,99,119,141 \
        --out-prefix ${sample} \
        -o ${outputdir}/Assembly/${sample}/megahit_out
    else
        echo "error! chose one for mode from PE or SE"
    fi

    if [[ -s ${outputdir}/Assembly/${sample}/megahit_out/${sample}.contigs.fa ]]
    then
        echo "${sample} megahit was done!" >> ${outputdir}/work_dir/metagenome_pipeline.log
    else
        echo "${sample} megahit error!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        exit
    fi

    MegahitendTime=`date +'%Y-%m-%d %H:%M:%S'`
    MegahitendTime_s=$(date --date="$MegahitendTime" +%s)
    echo -e "${sample}\tMegahit\t${threads}\t$MegahitstartTime\t$MegahitendTime\t"$((MegahitendTime_s-MegahitstartTime_s)) >> ${outputdir}/work_dir/metagenome_pipeline_time.log

}

function MetaWrap(){

    assembly=${outputdir}/Assembly/${sample}/megahit_out/${sample}.contigs.fa
    suffix=megahit


    if [ -d ${outputdir}/Binning/metawrap/${sample}_${suffix} ]; then
        rm -rf ${outputdir}/Binning/metawrap/${sample}_${suffix}
    fi
    mkdir -p ${outputdir}/Binning/metawrap/${sample}_${suffix}

    MetawrapstartTime=`date +'%Y-%m-%d %H:%M:%S'`
    MetawrapstartTime_s=$(date --date="$MetawrapstartTime" +%s)


    if [ $mode == "PE" ]; then
        metawrap binning \
        -o ${outputdir}/Binning/metawrap/${sample}_${suffix} \
        -t ${threads} \
        -a ${assembly} \
        --metabat2 \
        --maxbin2 \
        --concoct \
        ${outputdir}/QC/${sample}/${sample}_final_1.fastq \
        ${outputdir}/QC/${sample}/${sample}_final_2.fastq
    elif [ $mode == "SE" ]; then
        metawrap binning \
        -o ${outputdir}/Binning/metawrap/${sample}_${suffix} \
        -t ${threads} \
        -a ${assembly} \
        --metabat2 \
        --maxbin2 \
        --concoct \
        --single-end \
        ${outputdir}/QC/${sample}/${sample}_final.fastq
    else
        echo "error! chose one for mode from PE or SE"
    fi


    if [ -d "${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins" -a -d "${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins" -a -d "${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins" ];then
        echo "${sample}_${suffix} Metawrap was done!" >> ${outputdir}/work_dir/metagenome_pipeline.log
    else
        echo "${sample}_${suffix} Metawrap was WRONG!" >> ${outputdir}/work_dir/metagenome_pipeline.log
    fi

    MetawrapendTime=`date +'%Y-%m-%d %H:%M:%S'`
    MetawrapendTime_s=$(date --date="$MetawrapendTime" +%s)
    echo -e "${sample}_${suffix}\tMetawrap\t${threads}\t$MetawrapstartTime\t$MetawrapendTime\t"$((MetawrapendTime_s-MetawrapstartTime_s)) >> ${outputdir}/work_dir/metagenome_pipeline_time.log

    bin_counter=0

    if [ "$(ls -A ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins)" ];then
        for file in `ls ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins`
        do
        newfile=`echo $file|sed "s/bin/${sample}_${suffix}_concoct_bin/g"`
        mv ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins/$file ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins/$newfile
        done
        concoct_dir=${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins
        concoct_num=0
        for F in ${concoct_dir}/*; do
                SIZE=$(stat -c%s "$F")
                if (( $SIZE > 50000)) && (( $SIZE < 20000000)); then 
                    concoct_num=$((concoct_num +1))
                fi
        done

        if [[ ${concoct_num} > 0 ]]; then
            bin_counter=`expr ${bin_counter} + 1`
        fi
    fi

    if [ "$(ls -A ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins)" ];then
        for file in `ls ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins`
        do
        newfile=`echo $file|sed "s/bin/${sample}_${suffix}_maxbin2_bin/g"`
        mv ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins/$file ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins/$newfile
        done

        maxbin2_dir=${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins
        maxbin2_num=0
        for F in ${maxbin2_dir}/*; do
                SIZE=$(stat -c%s "$F")
                if (( $SIZE > 50000)) && (( $SIZE < 20000000)); then 
                    maxbin2_num=$((maxbin2_num +1))
                fi
        done

        if [[ ${maxbin2_num} > 0 ]]; then
            bin_counter=`expr ${bin_counter} + 2`
        fi
    fi

    if [ "$(ls -A ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins)" ];then
        for file in `ls ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins`
        do
        newfile=`echo $file|sed "s/bin/${sample}_${suffix}_metabat2_bin/g"`
        mv ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins/$file ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins/$newfile
        done
        metabat2_dir=${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins
        metabat2_num=0
        for F in ${metabat2_dir}/*; do
                SIZE=$(stat -c%s "$F")
                if (( $SIZE > 50000)) && (( $SIZE < 20000000)); then 
                    metabat2_num=$((metabat2_num +1))
                fi
        done

        if [[ ${metabat2_num} > 0 ]]; then
            bin_counter=`expr ${bin_counter} + 4`
        fi
    fi

    if [ $bin_refinement ]; then
        if [ -d ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement ]; then
            rm -rf ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement
        fi
        mkdir -p ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement
        startTime=`date +'%Y-%m-%d %H:%M:%S'`
        startTime_s=$(date --date="$startTime" +%s)

        if [ ${bin_counter} = 0 ];then
            echo "${sample}_${suffix} metawrap_bin_refinement no bins for use!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        if [ ${bin_counter} = 1 ];then
            echo "${sample}_${suffix} metawrap_bin_refinement only concoct bins for use, no need to run refinement, please check!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        if [ ${bin_counter} = 2 ];then
            echo "${sample}_${suffix} metawrap_bin_refinement only maxbin2 bins for use, no need to run refinement, please check!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        if [ ${bin_counter} = 4 ];then
            echo "${sample}_${suffix} metawrap_bin_refinement only metabat2 bins for use, no need to run refinement, please check!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        if [ ${bin_counter} = 3 ];then
            metawrap bin_refinement \
            -t ${threads} \
            --keep-ambiguous \
            -c 50 \
            -x 10 \
            -o ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement \
            -A ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins \
            -B ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins
            echo "${sample}_${suffix} metawrap_bin_refinement maxbin2 concoct bins for use!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        if [ ${bin_counter} = 5 ];then
            metawrap bin_refinement \
            -t ${threads} \
            --keep-ambiguous \
            -c 50 \
            -x 10 \
            -o ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement \
            -A ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins \
            -B ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins
            echo "${sample}_${suffix} metawrap_bin_refinement maxbin2 concoct bins for use!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        if [ ${bin_counter} = 6 ];then
            metawrap bin_refinement \
            -t ${threads} \
            --keep-ambiguous \
            -c 50 \
            -x 10 \
            -o ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement \
            -A ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins \
            -B ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins
            echo "${sample}_${suffix} metawrap_bin_refinement maxbin2 metabat2 bins for use!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        if [ ${bin_counter} = 7 ];then
            metawrap bin_refinement \
            -t ${threads} \
            --keep-ambiguous \
            -c 50 \
            -x 10 \
            -o ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement \
            -A ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins \
            -B ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins \
            -C ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins
            echo "${sample}_${suffix} metawrap_bin_refinement 3 software bins for use!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        endTime=`date +'%Y-%m-%d %H:%M:%S'`
        endTime_s=$(date --date="$endTime" +%s)
        echo -e "${sample}_${suffix}\tmetawrap_bin_refinement\t${threads}\t$startTime\t$endTime\t"$((endTime_s-startTime_s)) >> ${outputdir}/work_dir/metagenome_pipeline_time.log

        for file in `ls ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement/metawrap_50_10_bins`
        do
        newfile=`echo $file|sed "s/bin/${sample}_${suffix}_bin_refinement/g"`
        mv ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement/metawrap_50_10_bins/$file ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement/metawrap_50_10_bins/$newfile
        done

    fi

}

function GetBbins_Binrefinement(){

    suffix=megahit
    assembly=${outputdir}/Assembly/${sample}/megahit_out/${sample}.contigs.fa


    software=metawrap
    workDir=${outputdir}/Binning/${software}/${sample}_${suffix}/bin_refinement


    if [ "$(ls -A ${workDir}/metawrap_50_10_bins)" ];then
        if [ -d ${workDir}/B_bins_mapping ]; then
            rm -rf ${workDir}/B_bins_mapping
        fi
        mkdir -p ${workDir}/B_bins_mapping

        GetBbinsstartTime=`date +'%Y-%m-%d %H:%M:%S'`
        GetBbinsstartTime_s=$(date --date="$GetBbinsstartTime" +%s)
        
        awk 'BEGIN{FS=OFS="\t"} {if($2 >=90 && $3 <= 5) {print}}' ${outputdir}/Binning/${software}/${sample}_${suffix}/bin_refinement/metawrap_50_10_bins.stats > ${workDir}/checkm_drep_B.tsv
        awk -v add='.fa' 'BEGIN{FS=OFS="\t"} {if($2 >=90 && $3 <= 5) {print $1 add}}' ${outputdir}/Binning/${software}/${sample}_${suffix}/bin_refinement/metawrap_50_10_bins.stats > ${workDir}/checkm_drep_B_name.txt
        sed -i s"|bin|${sample}_${suffix}_bin_refinement|"g ${workDir}/checkm_drep_B_name.txt
        for Bbin in `cat ${workDir}/checkm_drep_B_name.txt`;do
            grep '>' ${workDir}/metawrap_50_10_bins/${Bbin} >>  ${workDir}/B_bins_mapping/all_contigs_name.txt
        done
        sort ${workDir}/B_bins_mapping/all_contigs_name.txt | uniq > ${workDir}/B_bins_mapping/all_contigs_name_sorted_uniq.txt
        perl /pipeline/co-assembly/pick_contigs_megahit.pl ${workDir}/B_bins_mapping/all_contigs_name_sorted_uniq.txt ${assembly} >  ${workDir}/B_bins_mapping/${sample}_B_bins.fa

        if [ -s "${workDir}/B_bins_mapping/${sample}_B_bins.fa" ];then
                echo "${sample} GetBbins was done!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        else
                echo "${sample} GetBbins was WRONG!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        fi
        
        GetBbinsendTime=`date +'%Y-%m-%d %H:%M:%S'`
        GetBbinsendTime_s=$(date --date="$GetBbinsendTime" +%s)
        echo -e "${sample}\tGetBbins\t${threads}\t$GetBbinsstartTime\t$GetBbinsendTime\t"$((GetBbinsendTime_s-GetBbinsstartTime_s)) >> ${outputdir}/work_dir/metagenome_pipeline_time.log
    else
            echo "${sample} GetBbins was WRONG! No dereplicated genomes!" >> ${outputdir}/work_dir/metagenome_pipeline.log
    fi
}


function GetBbinsUnmapped_Binrefinement_bwa(){

    suffix=megahit

    software=metawrap
    workDir=${outputdir}/Binning/${software}/${sample}_${suffix}/bin_refinement

    if [ -d ${workDir}/B_bins_mapping ]; then
        mkdir -p ${workDir}/B_bins_mapping
    fi

    if [ -s "${workDir}/B_bins_mapping/${sample}_B_bins.fa" ]; then

        if [ ! -s "${workDir}/B_bins_mapping/${sample}_B_bins.fa.amb" ];then
        bwa index -a bwtsw ${workDir}/B_bins_mapping/${sample}_B_bins.fa -p ${workDir}/B_bins_mapping/${sample}_B_bins.fa
        fi

        if [ $mode == "PE" ]; then
            BWAstartTime=`date +'%Y-%m-%d %H:%M:%S'`
            BWAstartTime_s=$(date --date="$BMAstartTime" +%s)

            bwa mem \
            -t ${threads} \
            -M \
            ${workDir}/B_bins_mapping/${sample}_B_bins.fa \
            ${outputdir}/QC/${sample}/${sample}_final_1.fastq \
            ${outputdir}/QC/${sample}/${sample}_final_2.fastq \
            -o ${workDir}/B_bins_mapping/${sample}_bwa.sam

            samtools view -bf 12 ${workDir}/B_bins_mapping/${sample}_bwa.sam > ${workDir}/B_bins_mapping/${sample}_bwa_unmapped.sam
            
            picard SamToFastq \
            I=${workDir}/B_bins_mapping/${sample}_bwa_unmapped.sam \
            F=${workDir}/B_bins_mapping/${sample}_B_bins_unmapped.1.fq \
            F2=${workDir}/B_bins_mapping/${sample}_B_bins_unmapped.2.fq \
            FU=${workDir}/B_bins_mapping/${sample}_B_bins_unmapped_U.fq

            if [[ -s ${workDir}/B_bins_mapping/${sample}_B_bins_unmapped.1.fq ]] && [[ -s ${workDir}/B_bins_mapping/${sample}_B_bins_unmapped.2.fq ]]
            then
                echo "${sample} GetBbinsUnmapped_bwa was done!" >> ${outputdir}/work_dir/metagenome_pipeline.log
            else
                echo "${sample} GetBbinsUnmapped_bwa error!" >> ${outputdir}/work_dir/metagenome_pipeline.log
                exit
            fi
            
            BWAendTime=`date +'%Y-%m-%d %H:%M:%S'`
            BWAendTime_s=$(date --date="$BWAendTime" +%s)
            echo -e "${sample}\tGetBbinsUnmapped_bwa\t${threads}\t$BWAstartTime\t$BWAendTime\t"$((BWAendTime_s-BWAstartTime_s)) >> ${outputdir}/work_dir/metagenome_pipeline_time.log
        elif [ $mode == "SE" ]; then
            BWAstartTime=`date +'%Y-%m-%d %H:%M:%S'`
            BWAstartTime_s=$(date --date="$BMAstartTime" +%s)

            bwa mem \
            -t ${threads} \
            -M \
            ${workDir}/B_bins_mapping/${sample}_B_bins.fa \
            ${outputdir}/QC/${sample}/${sample}_final.fastq \
            -o ${workDir}/B_bins_mapping/${sample}_bwa.sam

            samtools view -bf 12 ${workDir}/B_bins_mapping/${sample}_bwa.sam > ${workDir}/B_bins_mapping/${sample}_bwa_unmapped.sam
            
            picard SamToFastq \
            I=${workDir}/B_bins_mapping/${sample}_bwa_unmapped.sam \
            F=${workDir}/B_bins_mapping/${sample}_B_bins_unmapped.fq \

            if [[ -s ${workDir}/B_bins_mapping/${sample}_B_bins_unmapped.fq ]]
            then
                echo "${sample} GetBbinsUnmapped_bwa was done!" >> ${outputdir}/work_dir/metagenome_pipeline.log
            else
                echo "${sample} GetBbinsUnmapped_bwa error!" >> ${outputdir}/work_dir/metagenome_pipeline.log
                exit
            fi
            
            BWAendTime=`date +'%Y-%m-%d %H:%M:%S'`
            BWAendTime_s=$(date --date="$BWAendTime" +%s)
            echo -e "${sample}\tGetBbinsUnmapped_bwa\t${threads}\t$BWAstartTime\t$BWAendTime\t"$((BWAendTime_s-BWAstartTime_s)) >> ${outputdir}/work_dir/metagenome_pipeline_time.log
        else
            echo "error! chose one for mode from PE or SE"
        fi

    else
        if [ $mode == "PE" ]; then
            ln -s ${outputdir}/QC/${sample}/${sample}_final_1.fastq ${workDir}/B_bins_mapping/${sample}_needgroup_unmapped.1.fq
            ln -s ${outputdir}/QC/${sample}/${sample}_final_2.fastq ${workDir}/B_bins_mapping/${sample}_needgroup_unmapped.2.fq
        elif [ $mode == "SE" ]; then
            ln -s ${outputdir}/QC/${sample}/${sample}_final.fastq ${workDir}/B_bins_mapping/${sample}_needgroup_unmapped.fq
        else
            echo "error! chose one for mode from PE or SE"
        fi
        echo "${sample} GetBbinsUnmapped no B bins link raw QC, mark as need groups!!" >> ${outputdir}/work_dir/metagenome_pipeline.log
        recordTime=`date +'%Y-%m-%d %H:%M:%S'`
        echo -e "${sample}\tGetBbinsUnmapped_bwa\t${threads}\t$recordTime\t$recordTime\t0" >> ${outputdir}/work_dir/metagenome_pipeline_time.log
    fi

}



for ((i=0;i<${#name[@]};i++))
do
    sample=${name[${i}]}
    startpoint=${startpoints[${i}]}
    if [ $mode == "PE" ]; then
        rawread1=${pe1[${i}]}
        rawread2=${pe2[${i}]}
    elif [ $mode == "SE" ]; then
        rawread=${se[${i}]}
    else
        echo "error! chose one for mode from PE or SE"
    fi

    signal=0

    for var in ${array[@]}
    do
           case $var in
            Megahit)
            if [ ${startpoint} == "Megahit" ];then
                signal=1
            fi
            if [ $signal == 1 ];then
                Megahit
            fi
            ;;
            MetaWrap)
            if [ ${startpoint} == "MetaWrap" ];then
                signal=1
            fi
            if [ $signal == 1 ];then
                MetaWrap
            fi
            ;;
            GetBbins_Binrefinement)
            if [ ${startpoint} == "GetBbins_Binrefinement" ];then
                signal=1
            fi
            if [ $signal == 1 ];then
                GetBbins_Binrefinement
            fi
            ;;
            GetBbinsUnmapped_Binrefinement_bwa)
            if [ ${startpoint} == "GetBbinsUnmapped_Binrefinement_bwa" ];then
                signal=1
            fi
            if [ $signal == 1 ];then
                GetBbinsUnmapped_Binrefinement_bwa
            fi
            ;;
            *)  echo 'wrong array number, please check!'
            ;;
           esac
    done
done





