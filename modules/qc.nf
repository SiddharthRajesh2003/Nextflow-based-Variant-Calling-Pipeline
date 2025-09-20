#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process Quality_Control {
    tag "Quality Control with NanoPlot"
    publishDir "${qc_dir}", mode: 'copy'
    
    input:
    path fastq
    val qc_dir
    
    output:
    path "QC/NanoStats.txt"
    path "QC/*.png"
    path "QC/*.html"
    path "QC/*.log"
    
    script:
    """
    NanoPlot --fastq ${fastq} --outdir QC
    """
}