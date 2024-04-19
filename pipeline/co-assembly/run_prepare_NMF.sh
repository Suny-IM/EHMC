

dir=$1
sample_list=$2

mkdir -p ${dir}/kraken_mpa_percentage
mkdir -p ${dir}/kraken_mpa_percentage_g
mkdir -p ${dir}/percentage_g_nmf



for sample in `cat ${sample_list}`;do
    python kreport2mpa.py -r ${dir}/kraken_report/${sample}.kreport2 -o ${dir}/kraken_mpa_percentage/${sample}.mpa.percentage.txt --percentages --no-intermediate-ranks
done


for sample in `cat ${sample_list}`;do
    perl /pipeline/co-assembly/filter_for_g.pl ${dir}/kraken_mpa_percentage/${sample}.mpa.percentage.txt > ${dir}/kraken_mpa_percentage_g/${sample}.mpa.g.percentage.txt
done


for sample in `cat ${sample_list}`;do
    sed -i "1i #\t${sample}" ${dir}/kraken_mpa_percentage_g/${sample}.mpa.g.percentage.txt
done

echo -e "python combine_mpa.py -i \\" >> combine_mpa.g.sh
for sample in `cat  ${sample_list}`;do
    echo -e "${dir}/kraken_mpa_percentage_g/${sample}.mpa.g.percentage.txt \\" >> combine_mpa.g.sh
done
echo -e "-o ${dir}/kraken_mpa_percentage_g/combine.g.map.txt" >> combine_mpa.g.sh

sh combine_mpa.g.sh

grep -v Homo ${dir}/kraken_mpa_percentage_g/combine.g.map.txt > ${dir}/kraken_mpa_percentage_g/combine.g.map.nohomo.txt

ln -s ${dir}/kraken_mpa_percentage_g/combine.g.map.nohomo.txt ${dir}/percentage_g_nmf


