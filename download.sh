#!/bin/bash

#SBATCH -J DNA-seq
#SBATCH -p general
#SBATCH -o sra.txt
#SBATCH -e sra.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sidrajes@iu.edu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --time=40:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

base=/N/project/Krolab/Siddharth/Personal/DNA-seq
cd $base

module load sra-toolkit
prefetch SRR34291219 --max-size 200G -p sra_file/

fasterq-dump sra_file/SRR34291219 -O $base/fastq