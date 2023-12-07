#!/bin/sh

# Activate conda environment with kallisto-bustools
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate scRNA

wd='/data/brooksma/scRNA-seq/06_200310'

# Index and transcript to gene file location
idx_dir='/data/brooksma/Index/Mouse/scRNA/scEiaD'
