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
        def sample_name = bam.baseName.replaceAll('_trimmed', '')
        """
        # Set JAVA memory options (use ~75GB of your 100GB allocation)
        export JAVA_OPTS="-Xmx75g -XX:ParallelGCThreads=4"
        
        # Create temp directory for GATK
        mkdir -p ./tmp
        
        gatk --java-options "-Xmx75g -XX:ParallelGCThreads=4" MarkDuplicates \\
            -I ${bam} \\
            -O ${sample_name}_marked_duplicates.bam \\
            -M ${sample_name}_marked_duplicates_metrics.txt \\
            --CREATE_INDEX true \\
            --VALIDATION_STRINGENCY LENIENT \\
            --TMP_DIR ./tmp \\
            --MAX_RECORDS_IN_RAM 5000000
        """
}