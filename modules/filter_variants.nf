#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process Filter_Variants {
    tag "Filter high-quality variants"
    publishDir "${params.vcf}/filtered", mode: 'copy'
    
    input:
    tuple path (vcf), path (vcf_index)
    
    output:
    tuple path ("*_filtered.vcf.gz"), path ("*_filtered.vcf.gz.tbi")
    
    script:
    def sample_name = vcf.baseName.replaceAll('\\.merge_output\\.vcf.*', '')
    """
    # Filter variants with quality >= 20 and depth >= 10
    bcftools filter \\
        -i 'QUAL>=20 && (INFO/DP>=10 || FORMAT/DP>=10)' \\
        -O z \\
        -o ${sample_name}_filtered.vcf.gz \\
        ${vcf}

    # Index the filtered VCF
    tabix -p vcf ${sample_name}_filtered.vcf.gz
    """
}