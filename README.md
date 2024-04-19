# The construction of metagenome-assembled genomes by a combinational strategy of single and co-assembly


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
| samtools    | 1.15          | https://github.com/samtools/samtools                                                       |
| picard      | 2.26.11       | https://github.com/broadinstitute/picard                                                   |
| dRep        | 3.2.2         | https://github.com/MrOlm/drep                                                              |
| R           | 4.1.3         | https://www.r-project.org/                                                                 |
| checkm      | 1.1.3         | https://github.com/Ecogenomics/CheckM                                                      |
| r-nmf       | 0.21.0        | https://cran.r-project.org/web/packages/NMF/index.html                                     |


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

Notice: For all script under the single_sample_assembly, it is recommand to take the parameters and commands for reference not use the script directly since the directory structure will be different.

##### Step1 QC

**script: pipeline/single_sample_assembly/qc.sh**

The script "qc.sh" contains the parameters and commands for using bbduk to trim raw data, using bowtie2 to remove potential contamination sequence, using fastqc to make a QC report for clean data.

##### Step2 Assembly


**script: pipeline/single_sample_assembly/
run_assembly.sh**

The script "run_assembly.sh" contains the parameters and commands for using spades for assembly and using megahit for assembly.


##### Step3 Binning

**script: pipeline/single_sample_assembly/
run_binning.sh**

The script "run_binning.sh" contains the parameters and commands for using metawrap for binning and refinement, and using quast to check bins' quality.


###  Co-assembly part

#### Preparation for co-assembly

##### Step1 Run kraken

Kraken results were used to generated the matrix used by NMF analysis.

**script: pipeline/co-assembly/run_kraken.sh**

After all kraken were done, create a directory for preparing the matrix, for example named 'preapre_kraken'.

Creating a subdirectory named "Kreken_report" under "prepare_kraken", and put all samples 's kreport2 files which you want to run NMF analysis in it. 

Making a file list all samples ' name you want to run NMF analysis, named this file for example, "sample_list", one sample name for one line.

Then, you may run the script "pipeline/co-assembly/run_prepare_NMF.sh" like below,

```
sh run_prepare_NMF.sh prepare_kraken sample_list
```

There are two parameters need for run_prepare_NMF.sh.
Parameter1: the directory with absolute path you created for preparing the matrix, here for example "preapre_kraken".
Parameter2: the file contain all samples' name you want to do NMF analysis, for example "sample_list"

After run_prepare_NMF.sh was done, you will get a
 file named "percentage_g_nmf/combine.g.map.nohomo.txt" under the directory "prepare_kraken", this is the file we used for NMF analysis.
 
 Notice,:
 1) we only keep the genus data from kraken result, and manually remove the result for "k__Eukaryota|k__Metazoa|p__Chordata|c__Dipnotetrapodomorpha|c__Mammalia|o__Primates|f__Hominidae|g__Homo".
 2) run_prepare_NMF.sh contain the script "pipeline/co-assembly/filter_for_g.pl".

##### Step2 Run NMF analysis and pick 30% samples from each cluster

```
Rscript pipeline/co-assembly/run_nmf.R
```

The script run_nmf.R will use the file "combine.g.map.nohomo.txt", and generate two pdf files named "consensusmap.pdf" and "nmf_rank_survey.pdf", and a file named "combine.g.normalize.nohomo.map.txt".

In script run_nmf.R, we set the chosing rank from 2 to 25, but in real operation, we also use 26 to 51 and 51 to 65.

Chosing a suitable rank number by "consensusmap.pdf" and "nmf_rank_survey.pdf" image.

![image](https://github.com/Suny-IM/EHMC/assets/166774491/4aacedad-d574-4903-ba4d-0d7945e9b961) ![image](https://github.com/Suny-IM/EHMC/assets/166774491/93dff78b-2797-416b-895c-733e9afa2cda)



For two picture showing as example here, we think rank 9 will be good enough, we chose the rank based on the number that cophenetic start to drop, also take consensusmap into consideration.

After chosing the rank number, running the script /pipeline/co-assembly/get_rank.R for getting the NMF cluster information.

Using the script get random to pick 30% sample num for each cluster.

```
Rscript get_rank.R 9
```

The script "/pipeline/co-assembly/get_rank.R" will use the file "combine.g.normalize.nohomo.map.txt" which was generated by run_nmf.R, and need one 
parameter for the rank you choose, here for example, set to 9.

After get_rank.R was done, you will get 2 png file named "nsNMF_rank_normalize.png" and "nsNMF_rank_normalize_small.png", 1 pdf file named "heatmap_cluster.pdf" and 1 csv file named "expr_group.csv".

The "expr_group.csv" contain two colum, first was the sample name, the second was the cluster name. 


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

    ![image](https://github.com/Suny-IM/EHMC/assets/166774491/45e99758-49b9-42e5-bf25-113fd775545d)


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

A bins means the genomes meet the requirement that completeness greater than or equal to 95 and contamination less than or equal to 5 and the quality score (quality socre = completeness - 5 x contamination) greater than or equal to 50 .


##### Step5 Removing the reads which mapped to the dRep combined Abins

For all picked samples from each NMF cluster (for resource limitation, we randomly picked 30% samples from each cluster), running bowtie2 to remove the reads which mapped to the combined A bins fasta which named as "drep_dir/A_bins_combine/all_A_contigs_combine.fa".

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

##### Step1 Run the iteration of assembly and binning

Create a directory for iteration co-assembly operation, here for example named "mix_dir".

For each cluster of NMF result, create a subdirectory under the "mix_dir", such like "mix_dir/Cluster1", "mix_dir/Cluster2" and so on.

For each cluster, prepare a file named "sample_list" under the cluster directory, which should contain four columns, **seperated by tab**, 

column1_cluster: name of this cluster

column2_read1 : this cluster's read1 combined fastq1 file name with absolute path

column3_read2 : this cluster's read2 combined fastq2 file name with absolute path

column4_analysis_start_point: Megahit

Take "Cluster1" for example, the "mix_dir/Cluster1/sample_list" will be like,

>Cluster1 /path_to/Cluster1_combine_unmapped.1.fq Cluster1_combine_unmapped.2.fq Megahit

 
Prepare a iteration co-assembly running shell for each cluster, , take "Cluster1" for example,

```
sh /pipeline/co-assembly/metagenome_iteration_mix.sh \
mix_dir/Cluster1 \
30 \
```
There are two parameters are required for metagenome_iteration_mix.sh, please put them in stable order like the example shows.

Parameter1 : The operation dir for this cluster's iteration co assemly, for example "mix_dir/Cluster1".
Parameter2 : The threads used in analysis, here for example, set to 30, change it according to the resource available.

Notice:
metagenome_iteration_mix.sh  contain the script "/pipeline/co-assembly/assembly_binning.sh".
assembly_binning.sh contain the script "/pipeline/co-assembly/pick_contigs_megahit.pl".

metagenome_iteration_mix.sh will do iteration assembly and binning, until the iteration's number reach to 3 or no B bins were generated in this iteration.

B bins means the genomes meet the requirement that completeness greater than or equal to 90 and contamination less than or equal to 5.

If the last iteration number for Cluster1 was 3, there will be 3 directory named as "I1", "I2" and "I3" generated under "mix_dir/Cluster1/Iteration".

##### Step2 Run the dRep for each cluster's genomes

After each cluster's iteration assembly and binning is done, run dRep for each cluster itself, take Cluster1 for example, assumed the iteration number for Cluster1 is 3.

```
sh mix_drep_50_5_50.sh \
mix_dir/Cluster1 \
3 \
Cluster1 \
30
```
There are four parameters are required for mix_drep_50_5_50.sh, please put them in stable order like the example shows.

Parameter1 : The operation dir for this cluster's iteration co assemly, for example "mix_dir/Cluster1".
Parameter2 :  the last iteration number of this cluster, for example "3"
Parameter3 : this cluster's name
Parameter4 : The threads used in analysis, here for example, set to 30, change it according to the resource available.

After running mix_drep_50_5_50.sh, take Cluster1 and the iteration number is 3 for example, you will get "mix_dir/Cluster1/result_3/allbins/drep/dereplicated_genomes", the fasta files under this directory are the co-assembly genomes result for this cluster.

### The combination of single sample assembly's and co-assembly's genomes

#### Step1 Running dRep for the combination of single sample assembly's and co-assembly's genomes

Create a directory named "drep_mix_single" for running the dRep fot the combination of single sample assembly's genomes and co-assembly's genome.

Put all  species-level representative metagenome-assembled genomes of single sample assembly and all co-assembly's clusters' genomes under the subdirectory "drep_mix_single/bins".

Making a genomeInfo.csv file as  described above.

Running the dRep command, 
```
dRep dereplicate \
drep_mix_single \
-g drep_mix_single/bins/*fa \
--genomeInfo drep_mix_single/genomeInfo.csv \
-sa 0.95 \
-nc 0.30 \
-ms 10000 \
-strW 0 \
-centW 0 \
-p 30 \
-comp 50 \
-con 5
```

#### Step2 Finding the co-assembly improved genomes and new genomes compared to single-assembly

Running /pipeline/combination_of_single_sample_assembly_and_co-assembly/pick_mix_new_improve.sh under the directory "drep_mix_single/data_tables".

```
sh /pipeline/combination_of_single_sample_assembly_and_co-assembly/pick_mix_new_improve.sh
```

Notice:
1. /pipeline/combination_of_single_sample_assembly_and_co-assembly/pick_mix_new_improve.sh contain a lot of small scripts under the "/pipeline/combination_of_single_sample_assembly_and_co-assembly", before you run the command, make sure all these small scirpts are available in your system environment.
2. pick_mix_new_improve.sh will only work if your co-assembly genomes files' name contain "Cluster" and your single sample assembly genomes files' name do not contain "Cluster".

After the completion of "pick_mix_new_improve.sh", you will obtain three files containing genomes' names: bin_new_num1.txt, bin_improve_num2.txt, and bin_all_num3.txt. Here, num1, num2, and num3 represent actual numbers. Specifically, num1 refers to new co-assembly genomes compared to species-level representative metagenome-assembled genomes of single sample assembly, num2 refers to improved co-assembly genomes compared to species-level representative metagenome-assembled genomes of single sample assembly, and num3 refers to all retained co-assembly genomes.


### The annotation for genomes

Same annotation 
parameters are used for both single sample assembly genomes and co-assembly genomes, you can find the details in the scirpt "the_annotation_for_genomes/run_annotation.sh" and "the_annotation_for_genomes/run_gtdbtk.sh".

It is recommand to take the parameters and commands for reference not use the script directly since the directory structure will be different. 

## Statistical analysis and visualization
Statistical analysis and visualization were handled by scripting with R languages. These scripts were placed in "Scripts" directory. All related input data for statistical analysis and visualization are in "Pre-processed_Files" directory.


## KEGG_Pathway_Module_Info
The file "KEGG_Pathway_Module_Info.txt" used in analysis were placed in "KEGG_Pathway_Module_Info" directory.
