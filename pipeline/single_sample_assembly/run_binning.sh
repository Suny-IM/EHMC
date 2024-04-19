

##Run binning for pair end reads
metawrap binning \
-o ${outputdir}/Binning/metawrap/${sample}_${suffix} \
-t ${threads} \
-a ${assembly} \
--metabat2 \
--maxbin2 \
--concoct \
${outputdir}/QC/${sample}/${sample}_final_1.fastq \
${outputdir}/QC/${sample}/${sample}_final_2.fastq


##Run binning for single end reads
metawrap binning \
-o ${outputdir}/Binning/metawrap/${sample}_${suffix} \
-t ${threads} \
-a ${assembly} \
--metabat2 \
--maxbin2 \
--concoct \
--single-end \
${outputdir}/QC/${sample}/${sample}_final.fastq


##Run binning refinement
metawrap bin_refinement \
-t ${threads} \
--keep-ambiguous \
-c 50 \
-x 10 \
-o ${outputdir}/Binning/metawrap/${sample}_${suffix}/bin_refinement \
-A ${outputdir}/Binning/metawrap/${sample}_${suffix}/maxbin2_bins \
-B ${outputdir}/Binning/metawrap/${sample}_${suffix}/metabat2_bins \
-C ${outputdir}/Binning/metawrap/${sample}_${suffix}/concoct_bins


##Run quast for bins
${quast_path}/quast.py \
-t ${threads} \
-o ${outputdir}/Bin/${sample}/${binname}/quast  \
${binDir}/${bins}