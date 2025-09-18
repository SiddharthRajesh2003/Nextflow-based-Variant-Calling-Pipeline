#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process AnnotateVariants {
    tag "Annotate variants with VEP"
    publishDir "${params.vcf}/annotated", mode:'copy'

    input:
    tuple path(vcf), path(vcf_index)
    path vep_cache_dir
    path ref
    path ref_fai

    output:
    tuple path("*_annotated.vcf.gz"), path("*_annotated.vcf.gz.tbi"), emit: annotated_vcf
    path "*_vep_summary.html", emit: vep_report
    path "*_vep_stats.txt", emit: vep_stats
    path "*.log", emit: logs

    script:
    def sample_name = vcf.baseName.replaceAll(/(_filtered)?\.vcf(\.gz)?$/, '')
    def cache_dir = vep_cache_dir ? "--dir_cache ${vep_cache_dir}" : ""
    def fork = params.vep_threads ?: 8
    """
    vep --input_file ${vcf} \\
        --output_file ${sample_name}_annotated.vcf.gz \\
        --format vcf \\
        --vcf \\
        --compress_output gzip \\
        --stats_file ${sample_name}_vep_stats.txt \\
        --stats_html ${sample_name}_vep_summary.html \\
        --warning_file ${sample_name}_vep.log \\
        --fork ${fork} \\
        --species homo_sapiens \\
        --assembly ${params.genome_build ?: 'GRCh38'} \\
        ${cache_dir} \\
        ${ref} \\
        --offline \\
        --cache \\
        --merged \\
        --everything \\
        --canonical \\
        --biotype \\
        --check_existing \\
        --variant_class \\
        --regulatory \\
        --protein \\
        --symbol \\
        --numbers \\
        --domains \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --shift_hgvs 1 \\
        --verbose

    # Index the annotated VCF
    tabix -p vcf ${sample_name}_annotated.vcf.gz
    
    # Log completion
    echo "VEP annotation completed for ${sample_name}" > ${sample_name}_vep_completion.log
    """
}