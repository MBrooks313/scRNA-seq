#######################################
# This is an CellRanger analysis snakemake script.
# Written my Matthew J. Brooks on November 15th, 2022
# This runs on an HPC running SLURM
#######################################



#######################################
# Import config file and modules needed
#######################################

# Import modules
import glob
import os
import json
import pandas as pd

# Snakemake Base location
try:
	BASE_DIR=os.environ['BASE_DIR']
except KeyError:
	print("I can not locate your BASE_DIR directory.")
	pass

# Import configs
configfile: BASE_DIR + "/cr_config.json"

# wd = config['work_dir']
ref = config["cr_ref"]


########################################
# Import sample names
########################################

# Get samples base name with fastqs
fqs = glob.glob(BASE_DIR + "/**/*R[12]*fastq.gz", recursive=True)

fq_samples = [i.split("/")[-2] for i in fqs]

# Get unique sample names
samples = []

for x in fq_samples:
    if x not in samples:
        samples.append(x)


#############################################################
# List of directories needed and end point files for analysis
#############################################################

CR = expand("{sample}/outs/web_summary.html", sample=samples)


##############################
# Snakemake rules for analysis
##############################

localrules: all

rule all:
        input:  CR
        params:
                batch = config["job_all"]


rule cellranger:
   input:
            BASE_DIR + "/{sample}"
   output:
            "{sample}/outs/web_summary.html",
            "{sample}/outs/raw_feature_bc_matrix.h5",
            "{sample}/outs/filtered_feature_bc_matrix.h5"
   log:    "logs/cellranger.{sample}.log"
   version: config["cr_ver"]
   params:
           rulename = "cellranger",
           batch = config["job_cr"],
           ref = config["cr_ref"]
   shell: """
   module load cellranger/{version} || exit 1;
   cellranger count \
   --id {input} \
   --fastqs {input} \
   --transcriptome {params.ref} \
   --localcores $SLURM_CPUS_PER_TASK \
   --localmem 120 \
   --jobmode=slurm \
   --maxjobs=20
   """
