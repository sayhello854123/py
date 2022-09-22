#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=02:00:00
#PBS -l mem=64gb
#PBS -S /bin/bash
#PBS -N h5admirror
#PBS -o /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/make_hca_h5admirror.out
#PBS -e /beegfs/work/hd_wh241/MI_revisions/analysis/jobs/make_hca_h5admirror.err
#PBS -q short
#PBS -m bea
#PBS -M roramirezf@uni-heidelberg.de

source ~/.bashrc;
conda activate sc_analysis;
cd /beegfs/work/hd_wh241/MI_revisions;

$CONDA_PREFIX/bin/Rscript ./analysis/utils/h5ad_mirrors.R \
        --scell_data "/beegfs/work/hd_wh241/MI_revisions/ext_data/hca_seurat.rds" \
        --out_file "/beegfs/work/hd_wh241/MI_revisions/ext_data/hca_seurat.h5ad";