#!/bin/bash

###prepare
for bins in `ls ../dereplicated_genomes`;do
    grep -w ${bins} ../genomeInfo.csv >> drep_result.csv
done
grep Cluster drep_result.csv > drep_result_mix.csv
grep -v Cluster drep_result.csv > drep_result_single.csv

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/combine_drep_cluster_score.pl Sdb.csv  Cdb.csv > combine_S_C.csv

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/check_drep_cluster.pl combine_S_C.csv > pure_mix_cluster.txt

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/check_drep_cluster_all.pl combine_S_C.csv > all_mix_cluster.txt

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/check_drep_cluster_same_with_single.pl combine_S_C.csv > single_mix_cluster.txt

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/check_drep_cluster_single.pl combine_S_C.csv > pure_single_cluster.txt

num1=`less pure_mix_cluster.txt|wc -l`
mv pure_mix_cluster.txt pure_mix_cluster_${num1}.txt

for id in `cat pure_mix_cluster_${num1}.txt`;do
    grep -w ${id} combine_S_C.csv >> combine_S_C_${num1}cluster.csv
done

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/add_qualityfor_bins.pl genomeInformation.csv combine_S_C_${num1}cluster.csv > combine_S_C_${num1}cluster_add.csv

for id in `cat pure_mix_cluster_${num1}.txt`;do
    grep -w ${id} Wdb.csv >> Wdb_${num1}.csv
done

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/add_qualityfor_bins.pl genomeInformation.csv Wdb_${num1}.csv > Wdb_${num1}_add.csv

awk -F ',' '{print $1}' Wdb_${num1}.csv > mix_bin_${num1}.txt


num2=`less drep_result_mix.csv | wc -l`
awk -F ',' '{print $1}' drep_result_mix.csv > mix_bin_${num2}.txt

sort mix_bin_${num2}.txt | uniq > mix_bin_${num2}_sorted.txt
sort mix_bin_${num1}.txt | uniq > mix_bin_${num1}_sorted.txt

num3=`expr ${num2} - ${num1}`

diff mix_bin_${num2}_sorted.txt mix_bin_${num1}_sorted.txt | grep '<' > mix_bin_${num3}_sorted.txt

sed -i s'|< ||'g mix_bin_${num3}_sorted.txt

for bin in `cat mix_bin_${num3}_sorted.txt`;do
    grep -w ${bin} combine_S_C.csv >> combine_S_C_${num3}.csv
done

awk -F ',' '{print $2}' combine_S_C_${num3}.csv | sort | uniq > combine_S_C_${num3}_id

awk -F',' -v OFS=',' 'NR==FNR{a[$1]=$1}NR>FNR{if($2 in a){print $0}}' combine_S_C_${num3}_id combine_S_C.csv > combine_S_C_${num3}cluster.csv

for id in `cat combine_S_C_${num3}_id`;do
    grep -w ${id} Wdb.csv >> Wdb_${num3}.csv
done

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/add_qualityfor_bins.pl genomeInformation.csv Wdb_${num3}.csv > Wdb_${num3}_add.csv

grep -v Cluster combine_S_C_${num3}cluster.csv >> combine_S_C_${num3}cluster_single.csv

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/check_drep_improve_cluster_pickhigh.pl combine_S_C_${num3}cluster_single.csv > cluster_${num3}_single_score.txt

grep Cluster combine_S_C_${num3}cluster.csv > combine_S_C_${num3}cluster_mix.csv

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/check_drep_improve_cluster_pickmix.pl cluster_${num3}_single_score.txt  combine_S_C_${num3}cluster_mix.csv > combine_S_C_${num3}cluster_mix_higher.csv

perl /pipeline/combination_of_single_sample_assembly_and_co-assembly/add_qualityfor_bins.pl genomeInformation.csv combine_S_C_${num3}cluster_mix_higher.csv > combine_S_C_${num3}cluster_mix_higher_add.csv

num4=`less combine_S_C_${num1}cluster.csv | wc -l`

awk -F ',' '{print $1}' combine_S_C_${num1}cluster.csv >> bin_new_${num4}.txt

num5=`less combine_S_C_${num3}cluster_mix_higher.csv | wc -l`

awk -F ',' '{print $1}' combine_S_C_${num3}cluster_mix_higher.csv >> bin_improve_${num5}.txt

num6=`expr ${num4} + ${num5}`

cat bin_new_${num4}.txt bin_improve_${num5}.txt | sort | uniq > bin_all_${num6}.txt