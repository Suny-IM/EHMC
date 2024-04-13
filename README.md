# 【The construction of metagenome-assembled genomes by a combinational strategy of single and co-assembly 临时脚本名】


This directory contains scripts related to the manuscript "【文章名】".

## PREPARE

### Software information

| Software    | Version       | Availability                                                                               |
|-------------|---------------|--------------------------------------------------------------------------------------------|
| bbduk       | 39.01         | https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/ |
| bowtie2     | 2.3.5         | https://github.com/BenLangmead/bowtie2                                                     |
| fastqc      | 0.11.9        | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/                                 |
| kraken2     | 2.1.2         | https://ccb.jhu.edu/software/kraken2/                                                      |
| megahit     | 1.2.9         | https://github.com/voutcn/megahit                                                          |
| spades.py   | 3.13.0        | https://github.com/ablab/spades                                                            |
| quast.py    | 5.2.0         | https://github.com/ablab/quast                                                             |
| metawrap    | 1.3.2         | https://github.com/bxlab/metaWRAP                                                          |
| gtdbtk      | 2.1.1         | https://github.com/Ecogenomics/GTDBTk                                                      |
| prodigal    | 2.6.3         | https://github.com/hyattpd/Prodigal                                                        |
| tRNAscan-SE | 2.0.12        | https://github.com/UCSC-LoweLab/tRNAscan-SE                                                |
| rnammer     | 1.2           | https://services.healthtech.dtu.dk/services/RNAmmer-1.2/                                   |
| cmscan      | 1.1.4         | http://eddylab.org/infernal/                                                               |
| antismash   | 6.1.1         | https://github.com/antismash/antismash                                                     |
| diamond     | 2.0.14.152    | https://github.com/bbuchfink/diamond                                                       |
| pfam_scan   | 1.6           | http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/                                             |
| RGI         | 6.0.2         | https://github.com/arpcard/rgi                                                             |


Note: Suppose all needed software commands availabled in system environment variables, for these software commands will use directly without specify a software installed path in all the scripts in the pipeline directroy. Or, you may create a separate conda environment for each software, and activate the environment before you run the command and deactivate the environment at the end of the script.

### Database information

|        Database       |                    Version/Time stamp                    | Availability/Database website                                                         |
|:---------------------:|:--------------------------------------------------------:|---------------------------------------------------------------------------------------|
|          CARD         | prevalence-v4.0.0,   broadstreet-v3.2.5, ontology-v3.2.5 |                               https://card.mcmaster.ca/                               |
| card.json   (RGI_use) | 2023/1/28                                                |                          https://card.mcmaster.ca/latest/data                         |
|          VFDB         | 2022/11/11                                               |                               http://www.mgc.ac.cn/VFs/                               |
|          COG          |  2022/3                                                  |                       https://www.ncbi.nlm.nih.gov/research/cog/                      |
|          CAZy         | CAZyDB.08062022                                          |                                  http://www.cazy.org/                                 |
|          KEGG         | KEGG FTP   Release 2022-11-07                            |                              https://www.genome.jp/kegg/                              |
|          Pfam         | 2021/11/15                                               |                              http://pfam-legacy.xfam.org/                             |
|        BacMet2        | version 2.0                                              | http://bacmet.biomedicine.gu.se/download_temporary.html                               |
|       SwissProt       | SwissProt   2022/10/12                                   | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/ |
|       UniRef100       | Release:   2022_04, 12-Oct-2022                          |    https://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref100/    |
|        MetaCyc        | MetaCyc 25.1                                             |                                  https://metacyc.org/                                 |
|          Rfam         | Release 14.9                                             |                                   https://rfam.org/                                   |


## INSTRUCTION OF PIPELINE

There are four part of metagenomic analysis, including single sample assembly part, co-assembly part, the combination of single sample assembly's and co-assembly's result and the annotation for genomes.

### Single sample assembly part

#### The construction of metagenome-assembled genomes single sample assembly

##### Step1 QC

**script: QC.sh**

##### Step2 Assembly

**script: Run_assembly.sh**

##### Step3 Binning

**script: Run_binning.sh**

###  Co-assembly part

#### Preparation for co-assembly

##### Step1 Run kraken

Kraken result was used to generated the matrix used by NMF analysis

**script: run_kraken.sh**

Transfer report to mpa file, filter to keep genus data and filter out non human data. 

##### Step2 Run NMF analysis and pick 30% samples from each cluster

**script: Run_nmf.R**

Chosing a suitable rank by png, after chosing the rank, getting the cluster info of NMF result.
Seeting a seed value this time for helping to repeat the result.

Using the script get random to pick 30% sample num for each cluster.

##### Step3 Run dRep for single sample assembly's genome and get the combined A bins fasta of all the species-level representative metagenome-assembled genomes of single sample assembly part.

##### Prepare for dRep

1. Create a directory for operation the dRep of all genomes generated by single sample assembly parts which meets the requirement that completeness 
greater than or equal to 50 and contamination 
less than or equal to 5 and the quality score (quality socre = completeness - 5 x contamination) greater than or equal to 50 , for example, "path_to_drep".

2. Create a subdirectory named "bins", for example "path_to_drep/bins", which contain all the genomes fasta files generated by single sample assembly parts.

3. Preare a file named "genomeInfo.csv" under the "path_to_drep", for example "path_to_drep/genomeInfo.csv".

   genomeInfo.csv should contain three columns seperated by comma and contain the header which is "genome,completeness,contamination".

   You may get the completeness and contamination info from the single sample assembly's binning result.

   The corresponding genome info you have in "path_to_drep/bins" need to list in the file "path_to_drep/genomeInfo.csv".

   Example lines for genomeInfo.csv

   ![4440d2462214078ed64b26b1c621929d.png](en-resource://database/6991:1)

##### Run dRep

  **scirpt: pipeline/co-assembly/run_dRep.sh**

  **useage:**

```
sh run_drep.sh path_to_drep 30
```

The first parameter is the absolute path of the directory you create for operation of dRep.

The second parameter is the threads you need to run the dRep, take 30 for example, you may change it by the resource available.

After the success running of dRep, you may get a directory named "dereplicated_genome" under the "path_to_drep", which contains all the species-level representative metagenome-assembled genomes of single sample assembly part.


**script: /pipeline/co-assembly/get_Abins_after_dRep.sh**

**useage:** 

```
sh /pipeline/co-assembly/get_Abins_after_dRep.sh
```
Note : Running the script /pipeline/co-assembly/get_Abins_after_dRep.sh under the directory "path_to_drep".

This step will generated a A_bins_combine directory under the "path_to_drep", which may contain the combined A bins fasta of all species-level representative metagenome-assembled genomes of single sample assembly part, named as "drep_dir/A_bins_combine/all_A_contigs_combine.fa".

A bins means the genomes meet the requirement that completeness greater than or equal to 95 and contamination less than or equal to 5 and the quality score (quality socre = completeness - 5 x contamination).


##### Step5 Removing the reads which mapped to the dRep combined Abins

For all picked samples from each NMF cluster, running bowtie2 to remove the reads which mapped to the combined A bins fasta which named as "drep_dir/A_bins_combine/all_A_contigs_combine.fa".

For each sample, the bowtiw2 command will be like below, suppose the sample name is "sample1" and threads available are 20.
```
bowtie2 -p 20  \
-x drep_dir/A_bins_combine/all_A_contigs_combine.fa \
-1 dir/QC/sample1/sample1_final_1.fastq \
-2 dir/QC/sample1/sample1_final_2.fastq \
--un-conc dir/QC/sample1/remove_A/sample1_un_conc.fq \
--un dir/QC/sample/remove_A/sample1_un.fq
```

After running bowtie2, you will get two new fastq files under the "dir/QC/sample1/remove_A", take sample1 for example, that will be "sample1_un_conc1.fq" and "sample1_un_conc.2.fq". These new fastq files are what we need to pooling together. According to the cluster info of NMF analysis result, samples in the same cluster will be pooled together.


##### Step6 Pooling sample together by the cluster info of NMF result

Samples' fastq files are pooled together by the command `cat` , however, **please notice that you should cat the read1 fastq and read2 fastq in the same samples' order**.

For example, there are three samples (sample1, sample2, sample3) belongs to Cluster1 need to be pooled together, the command will be like below.

```
#read1 combination
cat sample1_un_conc.1.fq >> Cluster1_combine_unmapped.1.fq
cat sample2_un_conc.1.fq >> Cluster1_combine_unmapped.1.fq
cat sample3_un_conc.1.fq >> Cluster1_combine_unmapped.1.fq

#read2 combination
cat sample1_un_conc.2.fq >> Cluster1_combine_unmapped.2.fq
cat sample2_un_conc.2.fq >> Cluster1_combine_unmapped.2.fq
cat sample3_un_conc.2.fq >> Cluster1_combine_unmapped.2.fq
```

#### Run co-assembly

##### Run the iteration of assembly and binning

Create a directory for iteration co-assembly operation, here for example named "mix_dir".

For each cluster of NMF result, create a subdirectory under the "mix_dir", such like "mix_dir/Cluster1", "mix_dir/Cluster2" and so on.

For each cluster, prepare a file named "sample_list" under the cluster directory, which should contain four columns, seperated by tab,, 

>column1_cluster: name of this cluster
>column2_read1 : this cluster's read1 combined fastq1 file name with absolute path
>column3_read2 : this cluster's read2 combined fastq2 file name with absolute path
>column4_analysis_start_point: Megahit

Take "Cluster1" for example, the "mix_dir/Cluster1/sample_list" will be like,

>Cluster1	/path_to/Cluster1_combine_unmapped.1.fq	Cluster1_combine_unmapped.2.fq	Megahit
 


### The combination of single sample assembly's and co-assembly's genome


### The annotation for genomes

Same steps for both single sample assembly genomes and co-assembly genomes.


##### Run_annotation.sh




