

##Run kraken for paired end reads
kraken2 \
--db ${kraken_database} \
--report ${outputdir}/KrakenClean/${sample}/${sample}.kreport2 \
--threads ${threads} \
--paired ${outputdir}/QC/${sample}/${sample}_final_1.fastq \
${outputdir}/QC/${sample}/${sample}_final_2.fastq \
--use-names \
--output ${outputdir}/KrakenClean/${sample}/${sample}.kraken2

##Run kraken for single end reads
kraken2 \
--db ${kraken_database} \
--report ${outputdir}/KrakenClean/${sample}/${sample}.kreport2 \
--threads ${threads} \
--use-names \
--output ${outputdir}/KrakenClean/${sample}/${sample}.kraken2 \
${outputdir}/QC/${sample}/${sample}_final.fastq