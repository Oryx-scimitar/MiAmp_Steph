
configfile: "config.yaml"

rule main:
  input: expand("results/{sample}/{sample}.database.details.tsv", sample=config["samples"])

rule trim:
 output:
  read1 = "results/{sample}/{sample}.trimmed_read1.fastq",
  read2 = "results/{sample}/{sample}.trimmed_read2.fastq",
  singles = "results/{sample}/{sample}.trimmed_singles.fastq",
  log = "results/{sample}/{sample}.trimming_by_sickle.log"
 input:
  read1 =lambda wildcards: "fastq/{}_1.fastq.gz".format(config['reads'][wildcards.sample]),
  read2 =lambda wildcards: "fastq/{}_2.fastq.gz".format(config['reads'][wildcards.sample]),
 params:
  qual_threshold = config["read_quality_threshold"],
  len_threshold = config["read_trim_length_threshold"],
  folder = "results/{sample}"
 shell:
  r"""sickle pe -f {input.read1} -r {input.read2} -o {output.read1} -p {output.read2} -s {output.singles} -t sanger -q {params.qual_threshold} -l {params.len_threshold} > {output.log}
      fastqc -f fastq -o {params.folder} {input.read1}
      fastqc -f fastq -o {params.folder} {input.read2}
   """

rule overlap:
 output: "results/{sample}/{sample}.overlap_by_flash.log"
 input:
  read1 = "results/{sample}/{sample}.trimmed_read1.fastq",
  read2 = "results/{sample}/{sample}.trimmed_read2.fastq",
 params:
  min_overlap = config["minimum_overlap"],
  max_overlap = config["maximum_overlap"],
  directory = "results/{sample}",
  prefix = "{sample}"
 shell:
  r"""flash -m {params.min_overlap} -M {params.max_overlap} -O -o {params.prefix} -d {params.directory} --threads=1 {input.read1} {input.read2} > {output}
      fastqc -f fastq -o {params.directory} {params.directory}/{params.prefix}.extendedFrags.fastq
   """

rule sort:
 output: "results/{sample}/{sample}.clusters.stats.tsv"
 input: "results/{sample}/{sample}.overlap_by_flash.log"
 params:
  prefix = "{sample}",
  directory = "results/{sample}",
  pcr = lambda wildcards: config['pcr'][wildcards.sample],
  primer = lambda wildcards: config['primers'][config['pcr'][wildcards.sample]]
 shell:
  "perl scripts/sortPrimers.pl --sample={params.prefix} --work_dir={params.directory} --primer={params.pcr} --primer_seq={params.primer}"

rule blast:
 output: 
  blast1 = "results/{sample}/{sample}.clusters.database.blast",
 input: "results/{sample}/{sample}.clusters.stats.tsv"
 params:
  reference = lambda wildcards: config['database'][config['pcr'][wildcards.sample]],
  query = "results/{sample}/{sample}.clusters.filtered.fasta"
 shell:
  r"""
  blastn -db {params.reference} -query {params.query} -outfmt '6 qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen evalue bitscore' -out {output.blast1} -max_target_seqs 10
  """
   
   
rule analyse_blast:
 output: 
  database = "results/{sample}/{sample}.database.details.tsv",
 input: 
  blast1 = "results/{sample}/{sample}.clusters.database.blast",
 params:
  prefix = "{sample}",
  directory = "results/{sample}",
  taxonomy = config['taxonomy'],
 shell:
  r"""
  perl scripts/reportMapping.pl --sample={params.prefix} --taxonomy={params.taxonomy} --work_dir={params.directory} --blast_file={input.blast1} --output={output.database}
  """






  
