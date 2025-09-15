#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MarkDuplicates {
        tag "Marking duplicates with GATK"
        publishDir "${params.aligned}", mode: 'copy'

        input:
        tuple path(bam), path(bai)

        output:
        tuple path("*_marked_duplicates.bam"), path("*_marked_duplicates.bai")
        path "*_marked_duplicates_metrics.txt"

        script:
        def sample_name = bam.baseName.replaceAll('_trimmed.bam', '')
        """
        gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_name}_marked_duplicates.bam \
        -M ${sample_name}_marked_duplicates_metrics.txt \
        --CREATE_INDEX true
        """
}