#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.base="/N/project/Krolab/Siddharth/Personal/DNA-seq"
params.fastq="${params.base}/fastq/SRR34291219.fastq"
params.fastq_dir="${params.base}/fastq"
params.ref_dir="${params.base}/reference"
params.ref="${params.base}/reference/human_ref.fa"
params.results="${params.base}/results"
params.qc_before_trim="${params.results}/fastqc/raw"
params.qc_after_trim="${params.results}/fastq/trimmed"
params.aligned="${params.results}/bam"
params.vcf="${params.results}/vcf"
params.logs_dir="${params.base}/logs_dir"
params.apps="${params.base}/Apps"
params.model_dir="${params.apps}/clair3_models/r941_prom_sup_g5014"
params.platform = "ont"  // Change to "hifi" for PacBio
params.clair3_threads = 8
params.clair3_chunk_size = 5000000

// Additional useful parameters
params.skip_trimming = false
params.skip_alignment = false
params.skip_duplicate_marking = false
params.fallback_to_alignment = true  // Add this parameter
params.help = false
params.fallback_to_duplicate_marking = true  // Add this parameter for consistency

def helpMessage() {
    log.info"""
    ===============================================
     DNA-seq Analysis Pipeline
    ===============================================
    
    Usage:
      nextflow run main.nf [options]
    
    Required parameters:
      --fastq           Path to input FASTQ file
      --ref             Path to reference genome FASTA
    
    Optional parameters:
      --platform        Sequencing platform [ont/hifi] (default: ont)
      --clair3_threads  Threads for Clair3 (default: 16)
      --skip_trimming   Skip read trimming step
      --skip_alignment  Skip alignment (use existing BAMs)
      --skip_duplicate_marking  Skip duplicate marking step
      --fallback_to_duplicate_marking   Fall back to duplicate marking if marked BAMs not found (default: true)
      --results         Output directory (default: results)
      --fallback_to_alignment   Fall back to alignment if BAMs not found (default: true)
    
    """.stripIndent()
}

include { Index_Reference } from './modules/index_ref.nf'
include { Minimap_Index } from './modules/index_ref.nf'
include { Quality_Control as Quality_Control_Raw } from './modules/qc.nf'
include { Quality_Control as Quality_Control_Trimmed } from './modules/qc.nf'
include { TrimReads} from './modules/trim.nf'
include { AlignReads } from './modules/align.nf'
include { MarkDuplicates } from './modules/mkdp.nf'
include { VariantCalling } from './modules/variant_call.nf'
include { VCFStats } from './modules/call_stats.nf'
include { Filter_Variants } from './modules/filter_variants.nf'
include { MultiQC_Report } from './modules/multiqc.nf'
include { AnnotateVariants } from './modules/variant_annotate.nf'

// Function to check if we should actually skip alignment
def shouldSkipAlignment() {
    if (!params.skip_alignment) {
        return false
    }
    
    def bamDir = file(params.aligned)
    if (!bamDir.exists()) {
        if (params.fallback_to_alignment) {
            log.warn "BAM directory ${params.aligned} does not exist! Will perform alignment instead."
            return false
        } else {
            error "BAM directory ${params.aligned} does not exist and fallback_to_alignment is disabled!"
        }
    }
    
    // Check if there are actually BAM files
    def bamFiles = bamDir.listFiles().findAll { it.name.endsWith('.bam') }
    if (bamFiles.isEmpty()) {
        if (params.fallback_to_alignment) {
            log.warn "No BAM files found in ${params.aligned}! Will perform alignment instead."
            return false
        } else {
            error "No BAM files found in ${params.aligned} and fallback_to_alignment is disabled!"
        }
    }
    
    return true
}

// Function to check if we should actually skip duplicate marking
def shouldSkipDuplicateMarking() {
    if (!params.skip_duplicate_marking) {
        return false
    }
    
    def bamDir = file(params.aligned)
    if (!bamDir.exists()) {
        if (params.fallback_to_duplicate_marking) {
            log.warn "BAM directory ${params.aligned} does not exist! Will perform duplicate marking instead."
            return false
        } else {
            error "BAM directory ${params.aligned} does not exist and fallback_to_duplicate_marking is disabled!"
        }
    }
    
    // Check if there are actually duplicate-marked BAM files
    def dupMarkedFiles = bamDir.listFiles().findAll { 
        it.name.endsWith('.bam') && (it.name.contains('_marked') || it.name.contains('_dedup') || it.name.contains('duplicates_marked') || it.name.contains('marked_duplicates'))
    }
    if (dupMarkedFiles.isEmpty()) {
        if (params.fallback_to_duplicate_marking) {
            log.warn "No duplicate-marked BAM files found in ${params.aligned}! Will perform duplicate marking instead."
            return false
        } else {
            error "No duplicate-marked BAM files found in ${params.aligned} and fallback_to_duplicate_marking is disabled!"
        }
    }
    
    return true
}


// Workflow definition
workflow {
    if (params.help) {
    helpMessage()
    System.exit(0)
    }

    if (!params.fastq || !params.ref) {
    log.error "ERROR: Please provide required parameters --fastq and --ref"
    helpMessage()
    exit 1
    }

    // Input channels
    fastq_ch = Channel.fromPath(params.fastq)
    ref_ch   = Channel.fromPath(params.ref)
    model_dir_ch = Channel.fromPath(params.model_dir)

    // Quality control and reference indexing (always run)
    qc_ch = Quality_Control_Raw(fastq_ch, params.qc_before_trim)
    
    // Index_Reference now returns a tuple: (fai_file, dict_file)
    ref_idx_tuple = Index_Reference(ref_ch)
    ref_fai_ch = ref_idx_tuple.map { fai, _dict -> fai }  // Extract just the .fai file
    
    // Determine if we should actually skip alignment
    def actuallySkipAlignment = shouldSkipAlignment()
    
    // Create minimap index if we're doing alignment
    if (!actuallySkipAlignment) {
        mmi_idx = Minimap_Index(ref_ch)
    }

    // Handle trimming
    if (params.skip_trimming) {
        log.info "Skipping read trimming - using original FASTQ"
        trimmed_fastq_ch = fastq_ch
        // Skip trimmed QC if we skip trimming
        trimmed_qc_ch = Channel.empty()
    } else {
        trimmed_fastq_ch = TrimReads(fastq_ch)
        trimmed_qc_ch = Quality_Control_Trimmed(trimmed_fastq_ch, params.qc_after_trim)
    }

    // BAM handling - cleaner logic
    if (actuallySkipAlignment) {
        log.info "Using existing BAM files from: ${params.aligned}"
        
        // Create a channel for existing BAM files
        existing_bam_ch = Channel.fromPath("${params.aligned}/*.bam", checkIfExists: true)
            .map { bam ->
                log.info "Processing existing BAM file: ${bam.name}"
                def bai_paths = [
                    file("${bam}.bai"),                           // standard .bam.bai
                    file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")  // alternative .bai
                ]
                
                def bai = bai_paths.find { it.exists() }
                if (bai) {
                    log.info "Found index for ${bam.name}: ${bai.name}"
                    return tuple(bam, bai)
                } else {
                    log.error "Index file not found for ${bam.name}"
                    log.error "Looked for: ${bai_paths.collect{it.toString()}.join(', ')}"
                    error "Index file not found for ${bam.name}. Please ensure BAM files are properly indexed."
                }
            }
        
        aln_ch = existing_bam_ch
        
    } else {
        log.info "Performing alignment with minimap2"
        
        // Check for existing BAMs that we can reuse
        existing_bam_ch = Channel.empty()
        if (file(params.aligned).exists()) {
            existing_bam_ch = Channel.fromPath("${params.aligned}/*.bam", checkIfExists: false)
                .map { bam ->
                    def bai_paths = [
                        file("${bam}.bai"),
                        file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")
                    ]
                    
                    def bai = bai_paths.find { it.exists() }
                    if (bai) {
                        log.info "Found existing indexed BAM: ${bam.name}"
                        return tuple(bam, bai)
                    } else {
                        log.warn "Found BAM without index: ${bam.name} - will skip"
                        return null
                    }
                }
                .filter { it != null }
        }
        
        // Perform new alignment
        new_bam_ch = AlignReads(trimmed_fastq_ch, mmi_idx)
        
        // Use cached BAMs if available, otherwise use new BAMs
        aln_ch = existing_bam_ch.mix(new_bam_ch)
    }

    def actuallySkipDuplicateMarking = shouldSkipDuplicateMarking()

    // Handle duplicate marking - enhanced logic
    if (actuallySkipDuplicateMarking) {
        log.info "Using existing duplicate-marked BAM files from: ${params.aligned}"
        
        // Create a channel for existing duplicate-marked BAM files
        existing_dup_marked_ch = Channel.fromPath("${params.aligned}/*{_marked,_dedup,duplicates_marked}*.bam", checkIfExists: true)
            .map { bam ->
                log.info "Processing existing duplicate-marked BAM file: ${bam.name}"
                def bai_paths = [
                    file("${bam}.bai"),                           // standard .bam.bai
                    file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")  // alternative .bai
                ]
                
                def bai = bai_paths.find { it.exists() }
                if (bai) {
                    log.info "Found index for duplicate-marked ${bam.name}: ${bai.name}"
                    return tuple(bam, bai)
                } else {
                    log.error "Index file not found for duplicate-marked ${bam.name}"
                    log.error "Looked for: ${bai_paths.collect{it.toString()}.join(', ')}"
                    error "Index file not found for duplicate-marked ${bam.name}. Please ensure BAM files are properly indexed."
                }
            }
        
        md_ch = existing_dup_marked_ch
        
    } else {
        log.info "Performing duplicate marking"
        
        // Check for existing duplicate-marked BAMs that we can reuse in the same directory
        existing_dup_marked_ch = Channel.empty()
        if (file(params.aligned).exists()) {
            existing_dup_marked_ch = Channel.fromPath("${params.aligned}/*{_marked,_dedup,duplicates_marked}*.bam", checkIfExists: false)
                .map { bam ->
                    def bai_paths = [
                        file("${bam}.bai"),
                        file("${bam.toString().replaceAll(/\.bam$/, '.bai')}")
                    ]
                    
                    def bai = bai_paths.find { it.exists() }
                    if (bai) {
                        log.info "Found existing indexed duplicate-marked BAM: ${bam.name}"
                        return tuple(bam, bai)
                    } else {
                        log.warn "Found duplicate-marked BAM without index: ${bam.name} - will skip"
                        return null
                    }
                }
                .filter { it != null }
        }
        
        // Perform new duplicate marking
        new_dup_marked_ch = MarkDuplicates(aln_ch)
        md_output = new_dup_marked_ch  // Store the full output
        md_ch = md_output[0]  // Extract the BAM channel
        
        // Mix existing and new duplicate-marked BAMs if available
        if (existing_dup_marked_ch) {
            md_ch = existing_dup_marked_ch.mix(md_ch)
        }
    }

    // Variant calling - FIXED: Properly handle the reference index
    log.info "Starting variant calling with Clair3"
    vc_ch = VariantCalling(
        md_ch,
        ref_ch,
        ref_fai_ch,
        model_dir_ch
    )

    // Downstream analysis
    vcf_stats_ch = VCFStats(vc_ch[1])
    filtered_vcf_ch = Filter_Variants(vc_ch[1], vc_ch[2])

    annotated_vcf_ch = AnnotateVariants(filtered_vcf_ch,
        vep_cache_dir: params.vep_cache_dir ?: null,
        ref_ch,
        ref_fai_ch)

    // Apply collect to the channel variables, not directly on ProcessName.out
    ignore = qc_ch[2].collect()
    ignore = vcf_stats_ch[0].collect()
    ignore = trimmed_qc_ch[2].collect()
    ignore = annotated_vcf_ch.collect()

}