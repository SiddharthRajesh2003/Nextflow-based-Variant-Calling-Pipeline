#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process VCFStats {
    tag "VCF Statistics"
    publishDir "${params.vcf}", mode: 'copy'

    input:
    path vcf

    output:
    path "*_stats.txt"
    path "*_summary.txt"

    script:
    def sample_name = vcf.baseName.replaceAll('\\.merge_output\\.vcf.*', '')
    """
    bcftools stats ${vcf} > ${sample_name}_stats.txt

    echo "=== Variant Calling Summary for ${sample_name} ===" > ${sample_name}_summary.txt
    echo "Number of SNPs: \$(bcftools view -v snps ${vcf} | grep -v '^#' | wc -l)" >> ${sample_name}_summary.txt
    echo "Number of INDELs: \$(bcftools view -v indels ${vcf} | grep -v '^#' | wc -l)" >> ${sample_name}_summary.txt
    echo "Total variants: \$(bcftools view ${vcf} | grep -v '^#' | wc -l)" >> ${sample_name}_summary.txt
    echo "Generated: \$(date)" >> ${sample_name}_summary.txt
    """

}