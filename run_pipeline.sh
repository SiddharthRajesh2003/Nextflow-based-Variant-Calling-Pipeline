#!/bin/bash
#SBATCH -J DNA-seq
#SBATCH -p general
#SBATCH -o workflow_%j.txt                        # Fixed: Added %j for job ID
#SBATCH -e workflow_%j.err                        # Fixed: Added %j for job ID
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sidrajes@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1                           # Fixed: Changed from ntasks-per-node
#SBATCH --cpus-per-task=4                    # Added: CPUs for head job
#SBATCH --time=40:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

echo "Starting Nextflow pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"

# Set base directory
base=/N/project/Krolab/Siddharth/Personal/DNA-seq
cd $base

# Create output directory variable (fixed case)
OUTPUT_DIR=${base}/results                   # Fixed: Changed from lowercase 'output'
mkdir -p ${OUTPUT_DIR}/pipeline_info         # Create directory for reports

echo "Base directory: $base"
echo "Output directory: $OUTPUT_DIR"
export NXF_WORK="${base}/work"

# Load modules
module load fastqc
module load sra-toolkit
module load bcftools
module load minimap
module load samtools
module load gatk
module load java/17.0.7
module load conda

echo "Modules loaded successfully"

# Activate conda environment
echo "Activating conda environment..."
conda activate clair3

echo "Starting pipeline execution..."

# Run Nextflow pipeline (fixed duplicate report options)
nextflow run workflow.nf \
    -profile slurm \
    -resume \
    -with-timeline ${OUTPUT_DIR}/pipeline_info/execution_timeline_${SLURM_JOB_ID}.html \
    -with-report ${OUTPUT_DIR}/pipeline_info/execution_report_${SLURM_JOB_ID}.html \
    -with-trace ${OUTPUT_DIR}/pipeline_info/execution_trace_${SLURM_JOB_ID}.txt \
    -with-dag ${OUTPUT_DIR}/pipeline_info/pipeline_dag_${SLURM_JOB_ID}.svg

# Capture exit status
EXIT_STATUS=$?

echo ""
echo "Pipeline finished at $(date)"
echo "Exit status: $EXIT_STATUS"

if [ $EXIT_STATUS -eq 0 ]; then
    echo "SUCCESS: Pipeline completed successfully!"
    echo "Results available at: $OUTPUT_DIR"
    echo "Reports available at: ${OUTPUT_DIR}/pipeline_info/"
else
    echo "FAILED: Pipeline failed with exit status $EXIT_STATUS"
    echo "Check logs for details"
fi

exit $EXIT_STATUS