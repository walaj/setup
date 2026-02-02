#!/bin/sh
#SBATCH -p medium
#SBATCH -J mcmicro
#SBATCH -o mcmicro-%J.log
#SBATCH -t 0-72:00
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END     # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jeremiah_wala@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH -o hostname_mcmicro_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_mcmicro_%j.err                 # File to w

module purge
module load java/jdk-17.0.7
nextflow pull labsyspharm/mcmicro

## symlink
ln -s /home/jaw34/projects/jhu/params.yml params.yml
ln -s /home/jaw34/projects/jhu/markers.orion.csv markers.csv

## preconfig
/n/groups/lsp/mcmicro/tools/o2/config_post_reg.sh . > unmicst-memory.config

## submit
nextflow run labsyspharm/mcmicro --in . -profile O2,WSI,GPU --publish_dir_mode link -w /n/scratch/users/j/jaw34/work -c unmicst-memory.config -c /home/jaw34/projects/jhu/custom.config
