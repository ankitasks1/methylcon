# methylcon
Amplicon methylation data analysis wrapper script

The methylcon pipeline has been created using:

1. Python
2. Snakemake
3. Nextflow

For each version and language used a specific syntax and set of files required

##############################################

   # Python Version

##############################################


<ins>Run_syntax</ins>

<code>python methylcon.py --outdir out --softsheet /Users/ankitverma/Documents/Archivio2/Amplicon/methylcon/softsheet.csv --samplesheet /Users/ankitverma/Documents/Archivio2/Amplicon/methylcon/samplesheet.csv --genome /Users/ankitverma/Documents/Archivio2/Amplicon/methylcon/genome/ --genome_prefix hg38 --Rscript /usr/local/bin/Rscript --plotscript script.R</code>


<ins><b>cat samplesheet.csv</b></ins>

SRR11207817,SRR11207817_1.fastq.gz,SRR11207817_2.fastq.gz

SRR11207820,SRR11207820_1.fastq.gz,SRR11207820_2.fastq.gz


<ins><b>cat softsheet.csv</b></ins>

fastqc,/Users/ankitverma/Documents/Archivio2/Amplicon/methylcon/softwares/FastQC/
trimgalore,/Users/ankitverma/Documents/Archivio2/Amplicon/methylcon/softwares/TrimGalore-0.6.7/
cutadapt,/Users/ankitverma/miniconda3/bin/
bismark,/Users/ankitverma/Documents/Archivio2/Amplicon/methylcon/softwares/Bismark-0.22.3/
bowtie2,/Users/ankitverma/Documents/Archivio2/Amplicon/methylcon/softwares/bowtie2-2.5.0-macos-arm64/

genome folder will be created by default

##############################################

  # Snakemake Version

##############################################

<ins>Run syntax</ins>

<code>snakemake --cores 1</code>

or 

<code>snakemake -cores 1 -s methylcon.smk</code>

<ins><b>cat samplesheet.csv</b></ins>
N15_S7_L001,N15_S7_L001_R1_001.fastq.gz,N15_S7_L001_R2_001.fastq.gz
N18_S8_L001,N18_S8_L001_R1_001.fastq.gz,N18_S8_L001_R2_001.fastq.gz

Softwares can be given as full path or installed from conda/pip

and full path given can be changed in user case.

#Index must be generated before hand (ideally it should generate by default but there is a bug which need correctiion). No other bug is present and pipeline will run smoothly.
<code>snakemake -cores 1 -s methylcon.smk prepare_index</code>

#for Nextflow
pipeline generation in progress
