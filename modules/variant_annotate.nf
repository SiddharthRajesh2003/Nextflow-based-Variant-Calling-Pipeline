#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process AnnotateVariants {
    tag "Annotate variants with VEP"
    publishDir "${params.vcf}/annotated", mode:'copy'

    input:
    path vcf
    path vcf_index

    output:
    path "*.annotated.vcf.gz"
    path "*.html"
    path "*.annotated.vcf.gz.tbi"

    script:
    def sample_name = vcf.baseName.replaceAll('_filtered.vcf.gz', '')
    """
    vep --input_file ${vcf} \\
        --output_file ${sample_name}
    """
}