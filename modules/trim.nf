#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process TrimReads {
    tag "Reads trimming"
    publishDir "${params.fastq_dir}", mode:'copy'

    input:
    path fastq

    output:
    path "*_trimmed.fastq"

    script:
    """
    NanoFilt -q 15 -l 5000 --headcrop 50 --tailcrop 50 < ${fastq} > ${fastq.baseName}_trimmed.fastq
    """

}