#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process Index_Reference {
        tag "Indexing the reference genome"
        publishDir params.ref_dir, mode: 'copy'

        input:
        path ref

        output:
        path "*.fai"
        path "*.dict"

        script:
        """
        samtools faidx $ref
        gatk CreateSequenceDictionary -R $ref -O ${ref.baseName}.dict
        """
} 

// Minimap index for the reference genome

process Minimap_Index {
        tag "Indexing the reference genome for long read alignment"

        publishDir params.ref_dir, mode: 'copy'

        input:
        path ref

        output:
        path "*.mmi"

        script:
        """
        minimap2 -d ${ref.baseName}.mmi ${ref}
        """
}