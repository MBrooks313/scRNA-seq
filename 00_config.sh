#!/bin/sh

# Activate conda environment with kallisto-bustools installed
source <path_to_conda.sh_profile>
conda activate <conda_env>

# Working directory
wd=<absolute_path_to_working_directory>

# Index and transcript to gene file location
idx_dir=<absolute_path_to_index_directory>
