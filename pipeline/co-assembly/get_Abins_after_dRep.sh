dir=${PWD}
mkdir -p ${dir}/A_bins_combine/A_bins_new_name

for bins in `ls dereplicated_genomes`;do
    grep -w ${bins} genomeInfo.csv >> drep_result.csv
done

awk -F ',' '{if($2>=95 && $3<=5 && $2-5*$3>=50) print}'  drep_result.csv > drep_result_Abins.csv
awk -F ',' '{print $1}' drep_result_Abins.csv > A_bins_name

for bins in `cat A_bins_name`;do
    cp ${dir}/dereplicated_genomes/${bins} ${dir}/A_bins_combine/A_bins_new_name
done

for bins in `cat A_bins_name`;do
    sed -i s"|>|>${bins}|"g ${dir}/A_bins_combine/A_bins_new_name/${bins}
done

for bins in `cat A_bins_name`;do
    grep '>' ${dir}/A_bins_combine/A_bins_new_name/${bins} >>  ${dir}/A_bins_combine/all_A_contigs_name_new.txt
done

for bins in `cat A_bins_name`;do
    cat  ${dir}/A_bins_combine/A_bins_new_name/${bins} >>  ${dir}/A_bins_combine/all_A_contigs_combine.fa
done


bowtie2-build \
--threads 50 \
drep_dir/A_bins_combine/all_A_contigs_combine.fa \
drep_dir/A_bins_combine/all_A_contigs_combine.fa
