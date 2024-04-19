
##Run spades for paired end reads
spades.py \
--meta \
-1 ${outputdir}/QC/${sample}/${sample}_final_1.fastq \
-2 ${outputdir}/QC/${sample}/${sample}_final_2.fastq \
-t ${threads} \
-m 900 \
-o ${outputdir}/Assembly/${sample}/metaspades

##Run spades for single end reads
spades.py \
--meta \
-s ${outputdir}/QC/${sample}/${sample}_final.fastq \
-t ${threads} \
-m 900 \
-o ${outputdir}/Assembly/${sample}/metaspades


##Run megahit for paired end reads
megahit \
-1 ${outputdir}/QC/${sample}/${sample}_final_1.fastq \
-2 ${outputdir}/QC/${sample}/${sample}_final_2.fastq \
-t ${threads} \
--k-list 21,29,39,59,79,99,119,141 \
--out-prefix ${sample} \
-o ${outputdir}/Assembly/${sample}/megahit_out

##Run megahit for single end reads
megahit \
-r ${outputdir}/QC/${sample}/${sample}_final.fastq \
-t ${threads} \
--k-list 21,29,39,59,79,99,119,141 \
--out-prefix ${sample} \
-o ${outputdir}/Assembly/${sample}/megahit_out


