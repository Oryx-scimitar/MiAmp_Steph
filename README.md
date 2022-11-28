# MiAmp (Miseq Amplicon sequencing data analysis)

These snakemake pipeline and interative perl scripts is developed to analyse the targeted amplicon sequencing data (such as 16S/18S) using custom primers.


### Workflow
---

<img width="617" alt="Screenshot 2022-08-26 at 10 02 19" src="https://user-images.githubusercontent.com/8590103/186868476-4f504bb6-7352-4e6d-8960-76fc7ddb689f.png">


As shows in a workflow figure, it is split into two seperate analysis steps: Snakemake pipeline and creating summary tables

## 1. Snakemake pipeline (steps in green backgroud)
This pipeline analyse the raw sequencing data for each sample and produces the final results for each primer pairs used in the samples.

   - It starts with quality trimming of raw data. 
   - The high quality trimmed reads are then overlapped to generate extedened amplicon sequences.
   - Each extended sequence is then searched for PCR primers used. It searched for the forward and reverse primers, remove the excat primer sequences from reads to remove PCR sites. 
   
### Data received from MiSeq (FASTQ):

We receive data for each individual index and the folder names are given sample/library ID in samplesheet spreadsheet during sample submission. Each sample have two sequencing read files – forward and reverse. 

All sequencing facility has their own format of data files and folder structure. Data we receive from MiSeq through Edinburgh Genomics have two read files (fastq.gz) – forward and reverse and text file with the information of both reads (.count file)

  1.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_1.fastq.count*
  2.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_1.fastq.gz*
  3.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_2.fastq.count*
  4.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_2.fastq.gz*

If you want to check what these files have, simply use command less <file> command.
For example, `less 180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_1.fastq.count`

If you are using same index for multiple samples, we must make sure we assign samples with the correct files in the config file (see below).

Raw reads files must end with *_1.fastq.gz* and *_2.fastq.gz* as Snakefile recognises forward and reverse read by this extension OR you can change the extension of the input reads in Snakefile

### Conda environment:

We need following packages installed in environment: Perl, Blast, Flash, Sickle. 
There are two ways to setup conda environment. 1) using yml file 2) manully installing all packages. 
1. Using yml file: I have attached yml file. Following command will help to install from this file: `conda env create -f environment.yml`
2. Manually creating env and installing packages using following commands. “mhc” is the name of the environment, you can choice any name.
```
conda create --name amplicon
proceed ([y]/n)?
y
source activate amplicon
conda install -c bioconda snakemake
conda install -c conda-forge perl
conda install -c bioconda blast
conda install -c bioconda flash
conda install -c bioconda sickle-trim
conda install -c bioconda fastqc
```
Everytime when you want to run scripts or pipeline, use following command to activate environment first: 
`source activate amplicon`

  
### Folders and files: 
  
Create a main folder with project name and in this project folder should have following sub-folders:
- **fastq**: It should have all sample folders with paired end data
- **fasta**: 
    - It should have all primer sequences in fasta format for all primer pairs used and indexed database.
    - The available database fot he primers used in this study is the SILVA database. You can use your own database. 
    - To index database, use this command: `makeblastdb -in <database.fa> -dbtype nucl`
    - Primer sequences in fasta format. Make sure forward primers have “for” and reverse primers have “rev” in header. It doesn’t matter what else they are called in header but for and rev should be there. 
    - Make sure you remove all ambiguous IUPAC nucleotide and put all possible primer sequences in fasta file. Check out the primer sequences in fasta folder.
    - The primer sequences used in this study are available. 
- **results**: This directory will have all output files generated by pipeline
- **summary**: This folder will be used to create summary files for each sample
- **scripts**: All the perl scripts that will be used by Snakemake and other steps will be here
- **config.yaml**: It should be in your main project folder. Check out the attached example config file. You can modify it according to data you are analysing.
    - Make sure you don’t use sample names starts with numeric digits. Sample names must start with alphabets. If sample names start with digits, the simple trick is to put an alphabet in front of all samples. That will work.
- **Snakefile**: It should be in your main project folder.
- **samplesheet.txt**: This file is the list of samples.

The most important file to run pipeline is config.yaml where you are listing all samples. Make sure you are using right sample ids and locations of raw reads. 

### Running pipeline:
Before running pipeline, make sure you have following steps done.
  * In main folder, all subfolders are present.
  * In main folder, Snakefile and config file are present.
  * Config file must have correct information in correct format. Use my test config file to cross check the format. 
  * All reads should be in separate sample folders in fastq folder
  * Database and primer sequences should be in fasta format in fasta folder
  * Primer sequences should only have ATGC sequences and forward and reverse should be names as for and rev
  * Sequence database must be indexed. 
  * All the perl scripts must be in scripts folder

Once above steps are done, follow the commands below:

  1. To check if snakemake pipeline is ready to run and how many steps will be run, use this command:
         `snakemake -p -n`
	 
  2. To run pipeline
         `snakemake -p -j 3`
      
      Depending on the how many cores/processors you can use on system, change -j. 
      Alternatavely, you can submit individual jobs for each sample on cluster (Check out the script *submit_snakemake_eddie.pl*
  
  It will take a while to finish running pipeline depending on number of clusters you want to check for artefacts in config file and number of samples you are analysing. 
  
  Once the pipeline is finished running, there will be many files presence in result folder for each sample. These outputs are summarised into final tables in next phase of pipeline
 
## OUTPUT files:
There are multiple files generated during all steps. The details of each file is as below:
### Quality trimming:
* FASTQC reports are generated for raw reads. The file names are same as raw reads for each samples with *_fastqc.html* and *_fastqc.zip* extension. 
* Trimmed reads in fastq format. There are three trimmed fastq files created by sickle - *sample.trimmed_read1.fastq* and *sample.trimmed_read2.fastq* are paired end reads and *sample.trimmed_singles.fastq* is the file with only single end reads (it contains reads that passed filter in either the forward or reverse direction, but not the other).
* *sample.trimming_by_sickle.log* has information of number of trimmed reads. 

### Overlap paired end reads:
* *sample.extendedFrags.fastq* is the fastq of overlaped (extended) reads. 
* *sample.extendedFrags_fastqc.html* and *sample.extendedFrags_fastqc.zip* is the FASTQC report of extended reads. 
* *sample.notCombined_1.fastq* and *sample.notCombined_2.fastq* are paired end reads that failed to overlap each other and hence stored in seperate files.
* *sample.overlap_by_flash.log* has the information of number of reads overlapped etc. 
* *sample.hist, sample.hist.innie, sample.hist.outie, sample.histogram, sample.histogram.innie, sample.histogram.outie* are the histogram files that are generated by Flash. It has information of the length of overlapped reads. 
	
### Primer search and clustering
* *sample.primer.info.txt* has the information of primer sequences. Ideally, amplicon sequence should start and end with forward and reverse primers respectively. But due to size of amplicon (being shorter than sequencing read length) or PCR condition, sometimes the sequence doesn't start with forward primer nor ends with reverse primer. This file is the tab deliminated information of primer sequences and number of reads. 
* *sample.seq.len_his.txt* has the information of sequence length after trimming primers off the sequence. 
* *sample.clusters.fasta* has the sequences of all clusters. Clusters are the 100% identical sequences. The fasta header of each sequence represent the cluster details: unique id/read counts/read percentage/length of sequence/type of cluster. The types of cluster could be one of these:
	1. chimera: sequence that has detected as PCR chimera
	2. 1bpVariant: sequene which has 1 mismatch with higher abundant sequence
	3. ambiguous: sequence which has ambiguous basecall (N) 
	4. filters: sequence which is not identified as any of above. 
* *sample.clusters.details.tsv* has the details of each cluster sequence.
* *sample.clusters.stats.tsv* has the number of reads for each steps. The columns are sample, total reads with primers, total overlapped reads with primers, total single reads with primers, total clusters, total single copies, total chimeric reads, total reads with 1bp mismatch, total reads with low coverage, total reads with ambiguous calls, total filtered reads, total chimeric clusters, total clusters with 1bp mismatch, total clusters with low coverage, total clusters with ambiguous calls, total filtered clusters.
	
### Blast
* *sample.clusters.database.blast* is the output of blast on database sequence. The database is trimmed sesquences of SILVA database for the primers used in this study. The columns are qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen evalue bitscore. The details of the columns are found here: https://www.ncbi.nlm.nih.gov/books/NBK279684/
* *sample.clusters.nr.blast* is the output of blast on NCBI nr database. The columns are qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen evalue bitscore sscinames scomnames staxids stitle. The details of the columns are found here: https://www.ncbi.nlm.nih.gov/books/NBK279684/
* *sample.database.details.tsv* is the table of summary of database blast. 
* *sample.nr.details.tsv* is the table of summary of NCBI NR blast.

	
## 2. Create summary tables
   There are bunch of scripts developed to summarise the result of each sample and create final tables. This phage is more dependent on which PCR primers were used in terms of filtering out the results using read counts, length of amplicon and percent identity with database sequence. The script summaryTables.pl can create these tables in folder summary. 
   To run it:
   
   `perl scripts/summaryTables.pl example.samplesheet.txt`
   
summary.discarded.txt, summary.txt and summary_read_counts.tsv are the main tables that has the information of all samples in single file. While the summary folder has the summary tables and fasta file created for individual samples.

#Running remote blast on NCBI NR database:
   
   `perl scripts/remote.blast.pl example.samplesheet.txt`
   
   This script runs blast remotely on ncbi nr database for all filtered sequences of all samples. There will be nr.blast.tsv files created for individual sampels in summary folder. 
  
 

---
This pipeline is developed based on the following study:
The reference paper is: Deepali Vasoya, Andy Law, Paolo Motta, Mingyan Yu, Adrian Muwonge, Elizabeth Cook, Xiaoying Li, Karen Bryson, Amanda MacCallam, Tatjana Sitt, PhilipToye, Barend Bronsvoort, Mick Watson, W. Ivan Morrison and Timothy Connelley. **_"Rapid identification of bovine MHCI haplotypes in genetically divergent cattle populations Using Next-Generation Sequencing."_**
