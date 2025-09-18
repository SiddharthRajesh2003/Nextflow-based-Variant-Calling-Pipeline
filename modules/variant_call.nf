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
        path "${bam.baseName.replaceAll('_marked_duplicates', '')}_clair3/"
        path "${bam.baseName.replaceAll('_marked_duplicates', '')}_clair3/merge_output.vcf.gz"
        path "${bam.baseName.replaceAll('_marked_duplicates', '')}_clair3/merge_output.vcf.gz.tbi"
        
        script:
        def sample_name = bam.baseName.replaceAll('_marked_duplicates', '')
        """
        # Set locale to avoid locale errors
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8
        
        # Create output directory
        mkdir -p ${sample_name}_clair3

        # For female sample: exclude Y chromosome and problematic contigs
        # Create bed file for chromosomes 1-22, X (excluding Y, GL*, KI*, MT)

        awk '\$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X)\$/ {print \$1"\\t0\\t"\$2}' ${ref_fai} > female_chromosomes.bed
        
        # Debug: Show what regions we're processing
        
        echo "=== Chromosomes to process (female sample) ===" >&2
        cat female_chromosomes.bed >&2
        echo "Total regions: \$(wc -l < female_chromosomes.bed)" >&2
        
        # Run Clair3 with proper parameter format (using = signs)
        run_clair3.sh \\
            --bam_fn=${bam} \\
            --ref_fn=${ref} \\
            --threads=${params.clair3_threads} \\
            --platform=${params.platform} \\
            --model_path=${model_dir} \\
            --output=${sample_name}_clair3 \\
            --chunk_size=${params.clair3_chunk_size} \\
            --bed_fn=female_chromosomes.bed \\
            --haploid_precise \\
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