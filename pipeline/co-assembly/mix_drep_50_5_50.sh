#!/bin/bash
######################################################################
###This script is used for iteration metagenome drep
######################################################################


dir=$1
iteration_num=$2
groupID=$3
threads=$4

software="megahit"


dirnum=${iteration_num}

if [ -d ${dir}/result_${dirnum}/allbins/drep/bins ]; then
    rm -rf ${dir}/result_${dirnum}/allbins/drep/bins
fi
mkdir -p ${dir}/result_${dirnum}/allbins/drep/bins

echo "genome,completeness,contamination" > ${dir}/result_${dirnum}/allbins/drep/genomeInfo.csv


while [ ${iteration_num} -gt 1 ]
do
    cp ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/metawrap_50_10_bins.stats ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp
    sed -i s"|bin|I${iteration_num}_${software}_bin_refinement|"g ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp
    awk -v add='.fa' -v OFS=',' '{if($2>=50 && $3<=5 && $2-5*$3>=50) print $1 add,$2,$3}' ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp > ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/Allbins_info.csv
    awk -v add='.fa' '{if($2>=50 && $3<=5 && $2-5*$3>=50) print $1 add}' ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp > ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/Allbins_name.csv
    cat ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/Allbins_info.csv >> ${dir}/result_${dirnum}/allbins/drep/genomeInfo.csv
    for bins in `cat ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/Allbins_name.csv`;do
        ln -s ${dir}/Iteration/I${iteration_num}/Binning/metawrap/I${iteration_num}_${software}/bin_refinement/metawrap_50_10_bins/${bins} ${dir}/result_${dirnum}/allbins/drep/bins
    done
    iteration_num=`expr ${iteration_num} - 1`
done

if [ -d ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software} ]; then
    cp ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/metawrap_50_10_bins.stats ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp
    sed -i s"|bin|${groupID}_${software}_bin_refinement|"g ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp
    awk -v add='.fa' -v OFS=',' '{if($2>=50 && $3<=5 && $2-5*$3>=50) print $1 add,$2,$3}'  ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp > ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/Allbins_info.csv
    awk -v add='.fa' '{if($2>=50 && $3<=5 && $2-5*$3>=50) print $1 add}'  ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/metawrap_50_10_bins.stats.tmp > ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/Allbins_name.csv

    cat ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/Allbins_info.csv >> ${dir}/result_${dirnum}/allbins/drep/genomeInfo.csv

    for bins in `cat ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/Allbins_name.csv`;do
        ln -s ${dir}/Iteration/I${iteration_num}/Binning/metawrap/${groupID}_${software}/bin_refinement/metawrap_50_10_bins/${bins} ${dir}/result_${dirnum}/allbins/drep/bins
    done
else
    echo "${groupID} has problem!"
fi

cd ${dir}/result_${dirnum}/allbins/drep/bins/
rename I ${groupID}I *fa
cd -

sed -i s"|I|${groupID}I|"g ${dir}/result_${dirnum}/allbins/drep/genomeInfo.csv

dRep dereplicate \
${dir}/result_${dirnum}/allbins/drep \
-g ${dir}/result_${dirnum}/allbins/drep/bins/*fa \
--genomeInfo ${dir}/result_${dirnum}/allbins/drep/genomeInfo.csv \
-sa 0.95 \
-nc 0.30 \
-ms 10000 \
-strW 0 \
-centW 0 \
-p ${threads} \
-comp 50 \
-con 5


for bins in `ls ${dir}/result_${dirnum}/allbins/drep/dereplicated_genomes`;do
    grep -w ${bins} ${dir}/result_${dirnum}/allbins/drep/genomeInfo.csv >> ${dir}/result_${dirnum}/allbins/drep/drep_result.csv
done
