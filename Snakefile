SAMPLES = []

with open("samplesheet.csv") as myfile:
  myfile = myfile.read().strip().split('\n')
  for lines in myfile:
     lines = lines.split(',')
     SAMPLES.append(lines[0])

print(SAMPLES)

rule all:
 input:
   expand("data/{sample}_R1_001.fastq", sample=SAMPLES), 
   expand("data/{sample}_R1_001_fastqc.html", sample=SAMPLES), 
   expand("data/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
   expand("data/{sample}_R2_001.fastq", sample=SAMPLES),
   expand("data/{sample}_R2_001_fastqc.html", sample=SAMPLES),
   expand("data/{sample}_R2_001_fastqc.zip", sample=SAMPLES),
   expand("data/{sample}_R1_001_val_1.fq", sample=SAMPLES),
   expand("data/{sample}_R2_001_val_2.fq", sample=SAMPLES),
   expand("genome/Bisulfite_Genome/"),
   expand("data/{sample}_R1_001_val_1_bismark_bt2_pe.bam", sample=SAMPLES),
   expand("data/{sample}_bismark_bt2_pe.bam", sample=SAMPLES),
   expand("data/{sample}_bismark_bt2_pe.sort.bam", sample=SAMPLES),
   expand("data/{sample}_bismark_bt2_pe_sortedByReadname.bam", sample=SAMPLES),
   expand("data/{sample}_bismark_bt2_pe.sort.bam.bai", sample=SAMPLES),
   expand("data/{sample}_bismark_bt2_pe_sortedByReadname.bam.bai", sample=SAMPLES),

rule unzip:
 input:
  R1="data/{sample}_R1_001.fastq.gz",
  R2="data/{sample}_R2_001.fastq.gz"

 output:
  R1="data/{sample}_R1_001.fastq",
  R2="data/{sample}_R2_001.fastq"

 shell:
  """
  gunzip -c {input.R1} > {output.R1}
  gunzip -c {input.R2} > {output.R2}
  """

rule fastqc:
 input:
  R1="data/{sample}_R1_001.fastq",
  R2="data/{sample}_R2_001.fastq"

 output:
  html1="data/{sample}_R1_001_fastqc.html",
  zip1="data/{sample}_R1_001_fastqc.zip",
  html2="data/{sample}_R2_001_fastqc.html",
  zip2="data/{sample}_R2_001_fastqc.zip"

 conda:
  "fastqc"
 
 shell:
  """
  fastqc {input.R1} -o data/ 
  fastqc {input.R2} -o data/
  """


rule trimgalore:
 input:
  R1="data/{sample}_R1_001.fastq",
  R2="data/{sample}_R2_001.fastq"

 output:
  "data/{sample}_R1_001_val_1.fq",
  "data/{sample}_R2_001_val_2.fq"

 shell:
  """
  ~/Documents/Archivio2/Amplicon/methylcon/softwares/TrimGalore-0.6.7/trim_galore --phred33 --length 36 -q 20 --paired {input.R1} {input.R2} --path_to_cutadapt cutadapt -o data
  """

rule prepare_index:
 input:
  fastadir = "genome/"
 output:
  directory("genome/Bisulfite_Genome/")
 run:
  if os.path.exists(output[0]):
   pass
  else:
   shell("~/Documents/Archivio2/Amplicon/methylcon/softwares/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --verbose --path_to_aligner ~/Documents/Archivio2/Amplicon/methylcon/softwares/bowtie2-2.5.0-macos-arm64/ {input.fastadir}")

rule bismark_align:
 input:
  R1="data/{sample}_R1_001_val_1.fq",
  R2="data/{sample}_R2_001_val_2.fq"
 output:
  "data/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
 shell:
  """
  ~/Documents/Archivio2/Amplicon/methylcon/softwares/Bismark-0.22.3/bismark --score_min L,0,-0.6 --genome genome/ -1 {input.R1} -2 {input.R2} --bowtie2 --path_to_bowtie2 ~/Documents/Archivio2/Amplicon/methylcon/softwares/bowtie2-2.5.0-macos-arm64/ -o data/
  """

rule rename_bam:
 input:
  bam="data/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
 output:
  newbam="data/{sample}_bismark_bt2_pe.bam"
 shell:
  """
  cp {input.bam} {output.newbam}
  """

rule sortBAM:
 input:
  "data/{sample}_bismark_bt2_pe.bam"

 output:
  "data/{sample}_bismark_bt2_pe.sort.bam"

 shell:
  """
  samtools sort -o {output} {input}
  """

rule sortBAMbyname:
 input:
  "data/{sample}_bismark_bt2_pe.bam"

 output:
  "data/{sample}_bismark_bt2_pe_sortedByReadname.bam"

 shell:
  """
  samtools sort -n {input} -o {output}
  """

rule BAMindex:
 input:
  sortcoord="data/{sample}_bismark_bt2_pe.sort.bam",
  sortname="data/{sample}_bismark_bt2_pe_sortedByReadname.bam"

 output:
  incoord="data/{sample}_bismark_bt2_pe.sort.bam.bai",
  inname="data/{sample}_bismark_bt2_pe_sortedByReadname.bam.bai"

 shell:
  """
  samtools index {input.sortcoord} {output.incoord}
  samtools index {input.sortname} {output.inname}
  """

rule methylation_extract:
 input:
  name="data/{sample}_bismark_bt2_pe_sortedByReadname.bam"
 output:
  ""
 shell:
  """
  ~/Documents/Archivio2/Amplicon/methylcon/softwares/Bismark-0.22.3/bismark_methylation_extractor --bedGraph -p --include_overlap --zero_based --cutoff 10 {input.name} --samtools_path /usr/local/bin/samtools --cytosine_report --genome_folder genome/ -o data/
  """
