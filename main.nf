#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    =============================================================
     Metagenomics assembly pipeline v${workflow.manifest.version}
     Author: ${workflow.manifest.author}
    ============================================================
    Usage:
    The typical command for running the pipeline is as follows:

    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

if (params.csv != "") {
    // Parse CSV file of FASTQ files locations
    if (params.print_files) {
    Channel
      .fromPath(params.csv)
      .splitCsv(header:true, strip:true)
      .map{ row-> tuple(row.sampleId.strip(), file(row.forward), file(row.reverse)) }
      .println()
  }
  else {
    Channel
      .fromPath(params.csv)
      .splitCsv(header:true)
      .map{ row-> tuple(row.sampleId.strip(), file(row.forward), file(row.reverse)) }
      .set { input_fastq }
  }
}
else {
    error "Invalid input selected! Input must be converted CSV file."
}

// Run Kneaddata to remove low quality/human contaminated reads
process kneaddata {
  publishDir "$params.outdir/kneaddata/$sample_id", mode: 'copy', pattern: "fastqc/*.{html,zip}"
  publishDir "$params.outdir/kneaddata/$sample_id/log", mode: 'copy', pattern: "*.{log}"
  publishDir "$params.outdir/kneaddata/$sample_id/report", mode: 'copy', pattern: "*.{txt}"

  input:
    set val(sample_id), file(read_r1), file(read_r2) from input_fastq
    path genome from params.human_genome

  output:
    tuple val(sample_id), file("${sample_id}-clean_paired_{1,2}.fastq.gz"), file("${sample_id}-unpaired.fastq.gz") into filter_fq_assembly
    file("${sample_id}-clean.log") into fiter_log
    file("${sample_id}-kneaddata.txt") into kneaddata_log
    file("fastqc/*") into trim_fastqc
    //file("debug.txt") into kneaddata_debug

  script:
    """
    # Run kneaddata filtering tool
    kneaddata \
    --input ${read_r1} \
    --input ${read_r2} \
    --output . \
    --output-prefix ${sample_id}-clean \
    --sequencer-source TruSeq3 \
    --reference-db ${genome} \
    --threads 16 \
    --run-fastqc-start \
    --run-fastqc-end

    # Get readcount table
    kneaddata_read_count_table \
    --input ./ \
    --output ${sample_id}-kneaddata.txt

    # Compress mapped paired reads
    gzip ${sample_id}-clean_paired_1.fastq
    gzip ${sample_id}-clean_paired_2.fastq

    # Compress mapped unpaired reads
    gzip ${sample_id}-clean_unmatched_1.fastq
    gzip ${sample_id}-clean_unmatched_2.fastq

    # Concat for metaspades assembly
    cat ${sample_id}-clean_unmatched_1.fastq.gz ${sample_id}-clean_unmatched_2.fastq.gz > ${sample_id}-unpaired.fastq.gz
    """
}

// Functional profiling using assembly
process assembly {
  publishDir "$params.outdir/assembly/$sample_id/metaspades", pattern: "*-final.fasta.g", mode: 'copy'
  publishDir "$params.outdir/assembly/$sample_id/metaspades_quast", pattern: "metaspades_quast/", mode: 'copy'
  publishDir "$params.outdir/assembly/$sample_id/metaviralspades", pattern: "*-mvs-*.fasta.gz", mode: 'copy'
  publishDir "$params.outdir/assembly/$sample_id/metaviralspades_quast", pattern: "metaviralspades_quast/", mode: 'copy'
  publishDir "$params.outdir/assembly/$sample_id/", pattern: "*.{log,txt}", mode: 'copy'

  input:
    set sample_id, file(paired), file(unpaired) from filter_fq_assembly

  output:
    file("${sample_id}-spades.log") into metaspades_log
    //file("${sample_id}-contigs.fasta.gz") into metaspades_contigs
    file("${sample_id}-scaffolds.fasta.gz") into metaspades_scaffolds
    file("${sample_id}-contigs-final.fasta.gz") into metaspades_contigs_final
    file("metaspades_quast/*") into metaspades_quast
    file("${sample_id}-mvs-contigs.fasta.gz") into metaviralspades_contigs
    file("metaviralspades_quast/*") into metaviralspades_quast
    //file("debug.txt") into metaspades_debug

  script:
    """
    # Assembly with metaspades
    metaspades.py \
    -1 ${paired[0]} \
    -2 ${paired[1]} \
    -o metaspades/ \
    -s ${unpaired} \
    --threads 64

    # Rename contigs
    mv metaspades/contigs.fasta ${sample_id}-contigs.fasta
    mv metaspades/scaffolds.fasta ${sample_id}-scaffolds.fasta
    mv metaspades/spades.log ${sample_id}-spades.log

    # debugging
    #ls -lha > debug.txt

    # Remove short assemblies
    reformat.sh \
    in=${sample_id}-contigs.fasta \
    out=${sample_id}-contigs-final.fasta \
    minlength=1000

    # Get assembly quality information
    quast.py \
    ${sample_id}-contigs-final.fasta \
    -o metaspades_quast/

    # Assembly with metaviralspades
    metaviralspades.py \
    -1 ${paired[0]} \
    -2 ${paired[1]} \
    -o metaviralspades/ \
    -s ${unpaired} \
    --threads 64

    # Rename contigs
    mv metaviralspades/contigs.fasta ${sample_id}-mvs-contigs.fasta

    # Get assembly quality information
    quast.py \
    ${sample_id}-mvs-contigs.fasta \
    -o metaviralspades_quast/

    # Compress
    gzip ${sample_id}-contigs.fasta
    gzip ${sample_id}-scaffolds.fasta
    gzip ${sample_id}-contigs-final.fasta
    gzip ${sample_id}-mvs-contigs.fasta
    """
}
