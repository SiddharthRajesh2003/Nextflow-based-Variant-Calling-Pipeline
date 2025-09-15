#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process AlignReads {
        tag "Aligning the reads to the reference genome"

        publishDir "${params.aligned}", mode: 'copy'

        input:
        path fastq
        path mmi_index

        output:
        tuple path("*.bam"), path("*.bai")

        script:
        def sample_name = fastq.baseName.replaceAll('_trimmed.fastq', '')
        """
        minimap2 -ax map-ont -t 8 ${mmi_index} ${fastq} > ${sample_name}.sam
    
        samtools view -@ 8 -bS ${sample_name}.sam | \\
        samtools sort -@ 8 -o ${sample_name}.bam
        
        samtools index ${sample_name}.bam
        
        rm ${sample_name}.sam
        """
}