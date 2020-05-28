# MHC-ampliconSeq-analysis


These snakemake pipeline and interative perl scripts is developed to analyse the MHCI amplicon sequencing data or any targeted amplicon sequencing data (such as 16S/18S). 
The protocol used:
* PCR primers are designed to target specific exons in DNA or cDNA libraries. 
* Libraries amplified and sequenced using Miseq/Hiseq - <b>overlapping paired end reads</b>
* Steps includes
	* Quality trimming (using Sickle-trim and merging overlapping reads (Flash)
	* Identify primers sequences (Multiple primer pairs can be used)
	* Filter out PCR errors, PCR chimaeras, Low abundant sequences and unexpected length etc. 
	* Blast high confident sequences on database
	* Define new sequences. 
  * Create summary file
* Please check config file to check the parameters for snakemake and data

The reference paper is: Deepali Vasoya, Andy Law, Paolo Motta, Mingyan Yu, Adrian Muwonge, Elizabeth Cook, Xiaoying Li, Karen Bryson, Amanda MacCallam, Tatjana Sitt, PhilipToye, Barend Bronsvoort, Mick Watson, W. Ivan Morrison and Timothy Connelley. **_"Rapid identification of bovine MHCI haplotypes in genetically divergent cattle populations Using Next-Generation Sequencing."_**
