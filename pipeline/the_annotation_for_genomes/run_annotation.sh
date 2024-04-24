

##prodigal
prodigal \
-a ${outputdir}/Bin/${sample}/${bin}/prodigal_1/${bin}.faa \
-d ${outputdir}/Bin/${sample}/${bin}/prodigal_1/${bin}.fna \
-f gff \
-g 11 \
-p single \
-m \
-i ${binDir}/${bin}.fa \
-o ${outputdir}/Bin/${sample}/${bin}/prodigal_1/${bin}.gff

##CARD
diamond blastp \
-p ${threads} \
-d ${carddb} \
-q ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/CARD/${prodigal_num}/${bin}_card_diamond.txt \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40\
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

##VFDB
diamond blastp \
-p ${threads} \
-d ${vfdb} \
-q ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/VFDB/${prodigal_num}/${bin}_VFDB_diamond.txt \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40\
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

##COG
diamond blastp \
-p ${threads} \
-d ${cogdb} \
-q ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/COG/${prodigal_num}/${bin}_COG_diamond.txt \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40\
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

##CAZy
diamond blastp \
-p ${threads} \
-d ${cazydb} \
-q ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/CAZy/${prodigal_num}/${bin}_CAZy_diamond.txt \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40\
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp


##KEGG
diamond blastp \
-p ${threads} \
-d ${keggdb} \
-q ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/KEGG/${prodigal_num}/${bin}_KEGG_diamond.txt \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40 \
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

##SwissProt
diamond blastp \
-p ${threads} \
-d ${swissprotdb} \
-q ${outputdir}/Bin/${sample}/${bin}//${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/SwissProt/${prodigal_num}/${bin}_SwissProt_diamond.txt  \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40 \
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

##UniRef100
diamond blastp \
-p ${threads} \
-d ${uniref100db} \
-q ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/UniRef100/${prodigal_num}/${bin}_UniRef100_diamond.txt  \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40 \
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

##MetaCyc
diamond blastp \
-p ${threads} \
-d ${metacycdb} \
-q ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-o ${outputdir}/Bin/${sample}/${bin}/MetaCyc25_1/${prodigal_num}/${bin}_MetaCyc_diamond.txt  \
-e 1e-5 \
-k 1 \
--max-hsps 1 \
--id 40 \
--query-cover 40 \
--subject-cover 40 \
--outfmt 6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

##Pfam
pfam_scan.pl \
-fasta ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa \
-dir ${pfamdir} \
-outfile ${outputdir}/Bin/${sample}/${bin}/Pfam/${prodigal_num}/Pfam_diamond.txt \
-cpu ${thread}

##tRNAscan
tRNAscan_conf="/path_to_file/trnascanSE/bin/tRNAscan-SE.conf"
tRNAscan-SE \
-G \
-qQ \
-y \
-o ${outputdir}/Bin/${sample}/${bin}/tRNASCAN/${bin}.contigs.fa.final.txt \
-f ${outputdir}/Bin/${sample}/${bin}/tRNASCAN/${bin}.contigs.fa.tRNA.secondary.structures.txt \
-m ${outputdir}/Bin/${sample}/${bin}/tRNASCAN/${bin}.contigs.fa.statistics.summary \
-a ${outputdir}/Bin/${sample}/${bin}/tRNASCAN/${bin}.contigs.fa.tRNA.fasta \
-l ${outputdir}/Bin/${sample}/${bin}/tRNASCAN/${bin}.contigs.fa.log \
-c $tRNAscan_conf \
--thread ${threads} \
${workDir}/metawrap_50_10_bins/${bin}.fa


##Rfam
input=${workDir}/metawrap_50_10_bins/${bin}.fa
size=$(esl-seqstat $input | awk -F ':' '/Total.*residues.*/{print $2}' | sed 's/^[ \t]*//g')
num=$(echo "scale=4;$size*2.0/1000000" | bc)

${cmscan_path}/cmscan \
-Z $num --cut_ga --rfam --nohmmonly \
--tblout ${outputdir}/Bin/${sample}/${bin}/Rfam/${bin}.tblout \
--fmt 2 \
--cpu ${threads} \
--clanin ${Rfam_clanin_file} \
${Rfam_cm_file} $input > ${outputdir}/Bin/${sample}/${bin}/Rfam/${bin}.cmscan

awk 'BEGIN{OFS="\t";}{if(FNR==1) print "target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue"; if(FNR>2 && $20!="=" && $0!~/^#/) print $2,$3,$4,$10,$11,$12,$17,$18; }' ${outputdir}/Bin/${sample}/${bin}/Rfam/${bin}.tblout > ${outputdir}/Bin/${sample}/${bin}/Rfam/${bin}.tblout.final.xls


##Antismash
antismash \
--taxon bacteria \
-c ${threads} \
--cb-knownclusters \
--cb-general  \
--output-dir  ${outputdir}/Bin/${sample}/${bin}/Antismash/result \
--genefinding-tool prodigal \
${workDir}/metawrap_50_10_bins/${bin}.fa


##RNAmmer for bacteria
RNAmmer_S="bac"
RNAmmer_m="lsu,tsu,ssu"

cd ${outputdir}/Bin/${sample}/${bin}/RNAmmer &&
rnammer \
-S $RNAmmer_S \
-multi \
-m $RNAmmer_m \
-xml ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.xml \
-f ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.fasta \
-h ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.hmmreport \
-gff ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.gff2 \
${workDir}/metawrap_50_10_bins/${bin}.fa &&
cd -


##RNAmmer for archaea
RNAmmer_S="arc"
RNAmmer_m="lsu,tsu,ssu"

cd ${outputdir}/Bin/${sample}/${bin}/RNAmmer &&
rnammer \
-S $RNAmmer_S \
-multi \
-m $RNAmmer_m \
-xml ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.xml \
-f ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.fasta \
-h ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.hmmreport \
-gff ${outputdir}/Bin/${sample}/${bin}/RNAmmer/${bin}.RNAmmer.gff2 \
${workDir}/metawrap_50_10_bins/${bin}.fa &&
cd -


##Run RGI
awk '{if(/^>/ && NR==1){print $0} else if(/^>/ && NR>1){print "";print $0}else{printf $0}}' ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.faa > ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.tmp.faa
sed -i 's#*$##g' ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.tmp.faa

rgi main -i ${outputdir}/Bin/${sample}/${bin}/${prodigal_num}/${bin}.tmp.faa -t protein -o ${outputdir}/Bin/${sample}/${bin}/RGI/${prodigal_num}/${bin}_RGI --clean -n ${threads} --include_loose -a DIAMOND