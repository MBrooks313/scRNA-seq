#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 22:52:04 2022

@author: brooksma
"""

import glob
import pandas as pd
from os import path

work_dir = "/Volumes/data-1/scRNA-seq/06_200310"

##----------##
# Get a dataframe of sample directories and their fastq files

# Grab all fastq files 
fqs = glob.glob(work_dir + "/**/*R[12]*fastq.gz", recursive=True)

# Make a dict of the directories and fastqs
d = {'sample': [i.split("/")[-2] for i in fqs],
     'fastq': [i.split("/")[-1] for i in fqs]}

# Make a DataFrame from the dict
df = pd.DataFrame(d)

# Sort the DataFrame
df = df.sort_values(by = ['sample', 'fastq'])


##----------##
# Make the export analysis for each sample

samples = df['sample'].unique().tolist()

# Make an export doc for each sample
for sample in samples:
    
    # bash doc for sbatch
    file_nm = "01_" + sample + ".sh"
    
    with open(path.join(work_dir, "src", file_nm), "w") as file:
        
        # Common lines for analysis
        shebang = "#!/bin/bash"
        usage = "# sbatch -c 16 --mem=24g --time=8:00:00 " + file_nm
        config_file = "source 00_config.sh"
        cd_dir = "cd $wd"
        line1 = "kb count \\"
        line2 = "-i $idx_dir/idx.idx \\"
        line3 = "-g $idx_dir/t2g.txt \\"
        line4 = "-x 10XV3 \\"
        line5 = "-o kall_out/" + sample + " \\"
        line6 = "--h5ad \\"
        line7 = "--filter bustools \\"
        line8 = "--overwrite -t 16 \\"
        line9 = "--verbose \\"
        nl = ""  
        common_lines = [shebang, nl, usage, nl, config_file, nl, cd_dir, nl, 
                        line1, line2,line3, line4, line5, line6, line7, line8,
                        line9]
        
        # Get fastq files for the sample and join with directory
        tmp_fqs = df[df['sample'] == sample]['fastq'].tolist()
        tmp_fqs = [path.join(sample, i) for i in tmp_fqs]
        sample_fqs = ["%s \\" % i for i in tmp_fqs[:-1]]
        sample_fqs = sample_fqs + [tmp_fqs[-1]]
        
        # Join all lines and write to file
        lines = common_lines + sample_fqs
        file.writelines("%s\n" % l for l in lines)
        






    
    
