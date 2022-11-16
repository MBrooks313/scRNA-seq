#!/usr/bin/env python
# coding: utf-8


import glob
import os
import json
import pandas as pd


# Import configs
with open("cr_config.json", "r") as jsonfile:
    config = json.load(jsonfile)


# Get samples with fastqs
fqs = glob.glob(config['work_dir'] + "/**/*R[12]*fastq.gz", recursive=True)

fq_samples = [i.split("/")[-2] for i in fqs]


# Get unique samples names
samples = []

for x in fq_samples:
    if x not in samples:
        samples.append(x)


# Make an sbatch script for each sample
for sample in samples:
    
    # bash doc for sbatch
    file_nm = "01_cr_" + sample + ".sh"
    
    with open(os.path.join(config["work_dir"], "src", file_nm), "w") as file:
        
        # Common lines for analysis
        shebang = "#!/bin/bash"
        usage = "# sbatch --cpus-per-task=32 --mem=120G --time=24:00:00 --gres=lscratch:100 " + file_nm
        md = "mkdir -p " + config["work_dir"] + "/cellranger_out"
        cd = "cd " + config["work_dir"] + "/cellranger_out"
        mod = "module load cellranger/" + config['cr_ver']
        line1 = "cellranger count \\"
        line2 = "--id " + sample + " \\"
        line3 = "--transcriptome " + config['cr_ref'] + " \\"
        line4 = "--fastqs " + os.path.join(config["work_dir"], sample) + " \\"
        line5 = "--localcores $SLURM_CPUS_PER_TASK \\"
        line6 = "--localmem 120 \\"
        line7 = "--jobmode slurm \\"
        line8 = "--maxjobs 20"
        nl = ""  
        common_lines = [shebang, nl, 
                        usage, nl,
                        md, cd, nl,
                        mod, nl, 
                        line1, line2, line3, line4, line5, line6, line7, line8]
        
        # Join all lines and write to file
        lines = common_lines
        file.writelines("%s\n" % l for l in lines)

