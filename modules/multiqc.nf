#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MultiQC_Report {
    tag "MultiQC Report"
    publishDir "${params.results}/multiqc", mode: 'copy'
    
    input:
    path fastqc_files
    path stats_files
    
    output:
    path "multiqc_report.html"
    path "multiqc_data/"
    
    script:
    """
    multiqc . --title "Long-read Sequencing Analysis Report" \\
        --filename multiqc_report \\
        --force \\
        --verbose
    """
}