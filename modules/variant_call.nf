#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process VariantCalling {
        tag "Variant Calling with Clair3 for long reads"
        publishDir "${params.vcf}", mode: 'copy'
        
        input:
        tuple path(bam), path(bai)
        path ref
        path ref_fai
        path model_dir
        
        output:
        path "${bam.baseName}_clair3/"
        path "${bam.baseName}_clair3/merge_output.vcf.gz"
        path "${bam.baseName}_clair3/merge_output.vcf.gz.tbi"
        
        script:
        def sample_name = bam.baseName.replaceAll('_marked_duplicates', '')
        """
        # Set locale to avoid locale errors
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        # Create output directory
        mkdir -p ${sample_name}_clair3

        grep -v "^Y" ${ref_fai} | cut -f1,2 | awk 'BEGIN{OFS="\\t"} {print \$1, "0", \$2}' > exclude_y.bed
        
        # Run Clair3 with Python instead of PyPy
        run_clair3.sh \\
            -b ${bam} \\
            -f ${ref} \\
            --threads ${params.clair3_threads} \\
            --platform ${params.platform} \\
            --model_path ${model_dir} \\
            --output ${sample_name}_clair3 \\
            --include_all_ctgs \\
            --haploid_precise \\
            --bed_fn=exclude_y.bed \\
            --print_ref_calls \\
            --python=python3 \\
            --pypy=python3
        
        # Compress and index the VCF output
        if [ -f "${sample_name}_clair3/merge_output.vcf" ]; then
            bgzip -f ${sample_name}_clair3/merge_output.vcf
        fi
        
        if [ -f "${sample_name}_clair3/merge_output.vcf.gz" ]; then
            tabix -f -p vcf ${sample_name}_clair3/merge_output.vcf.gz
        fi
        """
}