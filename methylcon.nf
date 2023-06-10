params.input = "$projectDir/data/samplesheet.csv"
params.outdir = "$projectDir/allouts"
params.genomedir = "$projectDir/genome"
params.fasta = "ref.fa"
params.bismark_index = "$projectDir/genome/Bisulfite_Genome"


log.info """\
    M E T H Y L C O N - N F   P I P E L I N E
    =========================================
    Samplesheet Path        : ${params.input}
    Genome Index Directory  : ${params.genomedir}    
    Genome Fasta File       : ${params.fasta}
    Output Files Directory  : ${params.outdir}
    """
    .stripIndent()

process FASTQC {
    publishDir params.outdir, mode: 'copy'
    debug true
    input:
    tuple val(sampleId), path(read1), path(read2)

    output:
    path "fastqc_${sampleId}_logs"

    script:
    """
    mkdir fastqc_${sampleId}_logs
    echo "Quality Check: FastQC"
    echo $sampleId
    fastqc $read1 -o fastqc_${sampleId}_logs
    fastqc $read2 -o fastqc_${sampleId}_logs
    """
}

process TRIMGALORE {
    debug true

    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(sampleId), path(read1), path(read2)

    output:
    path "trimgalore_${sampleId}_logs"
    
    script:
    """
    echo "Trimming: TrimGalore/Cutadapt"
    /Users/ankitverma/Documents/tutorial/dollar_education/softwares/TrimGalore-0.6.7/trim_galore --phred33 --length 36 -q 20 --paired $read1 $read2 --path_to_cutadapt cutadapt -o trimgalore_${sampleId}_logs
    echo "Trimming done for $sampleId"
    """
}

process GENOME_INDEX {
    input:
    path "params.genomedir"

    //output is not required

    script:
    """
    echo "Generating Index: Reference Genome"
    echo "${params.fasta}"
    ~/Documents/tutorial/dollar_education/softwares/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --verbose --path_to_aligner ~/Documents/tutorial/dollar_education/softwares/bowtie2-2.5.0-macos-arm64/ ${params.genomedir}
    echo "Bismark Index generated"
    """
}

process BISMARK_ALIGN {
    debug true

    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(sampleId), path(read1), path(read2)
    path "trimgalore_${sampleId}_logs"

    output:
    path "bismark_${sampleId}_logs"    

    script:
    """
    echo "Trimmed Read 1: trimgalore_${sampleId}_logs/${sampleId}_1_val_1.fq.gz"
    echo "Trimmed Read 2: trimgalore_${sampleId}_logs/${sampleId}_2_val_2.fq.gz"
    echo "Alignment: Bismark"
    ~/Documents/tutorial/dollar_education/softwares/Bismark-0.22.3/bismark --score_min L,0,-0.6 --genome /Users/ankitverma/Documents/tutorial/nextflow/genome/ -1 trimgalore_${sampleId}_logs/${sampleId}_1_val_1.fq.gz -2 trimgalore_${sampleId}_logs/${sampleId}_2_val_2.fq.gz --bowtie2 --path_to_bowtie2 ~/Documents/tutorial/dollar_education/softwares/bowtie2-2.5.0-macos-arm64/ -o "bismark_${sampleId}_logs" 
    echo "Alignment done for $sampleId"
    """
}

process BAM_SORT {
    debug true

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(read1), path(read2)
    path "bismark_${sampleId}_logs"

    output:
    path "bismark_${sampleId}_logs"

    script:
    """
    echo "Renaming BAM ${sampleId}"
    cp bismark_${sampleId}_logs/${sampleId}_1_val_1_bismark_bt2_pe.bam bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.bam
    echo "Sorting BAM ${sampleId}"
    samtools sort -o bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.sorted.bam bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.bam
    echo "BAM sorted ${sampleId}"
    """
}

process BAM_SORT_READNAME {
    debug true

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(read1), path(read2)
    path "bismark_${sampleId}_logs"

    output:
    path "bismark_${sampleId}_logs"

    script:
    """
    echo "Sorting BAM by Readname ${sampleId}"
    samtools sort -n bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.bam -o bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.sortedByReadname.bam
    echo "BAM sorted by readname ${sampleId}"
    """
}

process BAM_INDEX {
    debug true

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(read1), path(read2)
    path "bismark_${sampleId}_logs"

    output:
    path "bismark_${sampleId}_logs"

    script:
    """
    echo "Indexing BAM ${sampleId}"
    samtools index bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.sorted.bam
    echo "Indexing ${sampleId}"
    """
}

process BAM_READNAME_INDEX {
    debug true

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sampleId), path(read1), path(read2)
    path "bismark_${sampleId}_logs"

    output:
    path "bismark_${sampleId}_logs"

    script:
    """
    echo "Indexing BAM ${sampleId}"
    samtools index bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.sortedByReadname.bam
    echo "Indexing ${sampleId}"
    """
}


process BISMARK_METHYLCALL {
    debug true

    publishDir params.outdir, mode: 'copy'
    input:
    tuple val(sampleId), path(read1), path(read2)
    path "bismark_${sampleId}_logs"

    output:
    path "bismark_${sampleId}_logs"

    script:
    """
    echo "Methylation calling: Bismark"
    ~/Documents/tutorial/dollar_education/softwares/Bismark-0.22.3/bismark_methylation_extractor --bedGraph -p --include_overlap --zero_based --cutoff 10 bismark_${sampleId}_logs/${sampleId}_bismark_bt2_pe.bam --samtools_path /usr/local/bin/samtools --cytosine_report --genome_folder /Users/ankitverma/Documents/tutorial/nextflow/genome/ -o bismark_${sampleId}_logs
    echo "Methylation calling done for $sampleId"
    """
}


process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path params.outdir

    output:
    path 'multiqc_report.html'

    script:
    """
    echo "Performing multiqc check"
    multiqc ${params.outdir} --verbose  --interactive --force
    """
}

workflow {
  Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sampleId, file(row.read1), file(row.read2)) }
        .set{ reads_ch }

   Channel
       .fromPath(params.fasta)
       .set{ fasta_ch }

   // Check if the output directory already exists
   def indexExists = file(params.bismark_index).exists()

   // Run GENOME_INDEX process only if the index doesn't exist
   if (!indexExists) {
       genome_index_ch = GENOME_INDEX(fasta_ch)
       println "Index were not found. So generating index"
       println "Finished Genome Index generation"
   }


  fastqc_ch = FASTQC(reads_ch)
  trim_ch = TRIMGALORE(reads_ch)
  //genome_index_ch = GENOME_INDEX(fasta_ch)
  bismark_align_ch = BISMARK_ALIGN(reads_ch, trim_ch)
  bam_sort_ch = BAM_SORT(reads_ch, bismark_align_ch)
  bam_sort_name_ch = BAM_SORT_READNAME(reads_ch, bam_sort_ch)
  bam_index_ch = BAM_INDEX(reads_ch, bam_sort_ch)
  bismark_methylcall_ch = BISMARK_METHYLCALL(reads_ch, bam_sort_name_ch)
  multiqc_ch = MULTIQC(bismark_methylcall_ch)
}

