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
params.clair3_threads = 1

// Additional useful parameters
params.skip_trimming = false
params.skip_alignment = false
params.skip_duplicate_marking = false
params.help = false

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
      --results         Output directory (default: results)
      --skip_duplicate_marking  Skip duplicate marking step
    
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



// Workflow definition
workflow {
    // Input channels
    fastq_ch = Channel.fromPath(params.fastq)
    ref_ch   = Channel.fromPath(params.ref)
    model_dir_ch = Channel.fromPath(params.model_dir)

    // Quality control and reference indexing (always run)
    qc_ch = Quality_Control_Raw(fastq_ch, params.qc_before_trim)
    ref_idx = Index_Reference(ref_ch)
    
    // Only create minimap index if we're doing alignment
    if (!params.skip_alignment) {
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

    // BAM handling - improved logic with better debugging
    if (params.skip_alignment) {
        log.info "Skipping alignment - looking for existing BAM files in: ${params.aligned}"
        
        // First, let's see what files actually exist
        def bamDir = file(params.aligned)
        if (!bamDir.exists()) {
            error "BAM directory ${params.aligned} does not exist!"
        }
        
        // Create a channel for existing BAM files with better error handling
        existing_bam_ch = Channel.fromPath("${params.aligned}/*.bam", checkIfExists: false)
            .ifEmpty { 
                log.error "No BAM files found in ${params.aligned}"
                log.error "Available files in directory:"
                try {
                    file(params.aligned).listFiles().each { 
                        log.error "  - ${it.name}" 
                    }
                } catch (Exception e) {
                    log.error "Could not list directory contents: ${e.message}"
                }
                error "No BAM files found in ${params.aligned} but alignment was skipped" 
            }
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

        new_bam_ch = AlignReads(trimmed_fastq_ch, mmi_idx)
        
        // Use cached BAMs if available, otherwise use new BAMs
        aln_ch = existing_bam_ch.mix(new_bam_ch)
    }

    // Mark duplicates (conditional)
    if (params.skip_duplicate_marking) {
        log.info "Skipping duplicate marking as requested"
        md_ch = aln_ch
    } else {
        log.info "Marking duplicates"
        md_ch = MarkDuplicates(aln_ch)
    }

    // Prepare reference files for variant calling
    ref_fai_ch = ref_idx[0]

    // Variant calling
    log.info "Starting variant calling with Clair3"
    vc_ch = VariantCalling(
        md_ch,
        ref_ch,
        ref_fai_ch,
        model_dir_ch
    )

    // Downstream analysis
    vcf_stats_ch = VCFStats(vc_ch[1])
    Filter_Variants(vc_ch[1], vc_ch[2])

    // Apply collect to the channel variables, not directly on ProcessName.out
    ignore  = qc_ch[2].collect()
    ignore = vcf_stats_ch[0].collect()
    ignore = trimmed_qc_ch[2].collect()
}