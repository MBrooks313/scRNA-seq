#!/bin/sh

# Run with: sbatch --time=48:00:00 cr_submit_snakemake.sh

#####################################
# This script is the CellRnager submit script for a snakemake pipeline.
# This pipeline was created by Matthew J Brooks in November 2022
# This pipeline adapted to run on HPCs running SLURM
# This requires the snakefile CellRanger.py and cr_config.json
#####################################

# Load module
module load python/3.7

# Export variables
NOW=$(date +"%Y%m%d")
# NOW='20210719'
export BASE_DIR="/data/brooksma/scRNA-seq/06_200310"
export WORK_DIR=${BASE_DIR}/cellranger_out/${NOW}
SNAKEFILE=${BASE_DIR}/CellRanger.py

# Make result directories and change into result directory
mkdir -p ${WORK_DIR}/logs
cd $WORK_DIR

# Snakemake command
echo "Get ready for snakemake..." >> logs/snakemake.%j.o
snakemake\
	--directory $WORK_DIR \
	--snakefile $SNAKEFILE \
	--jobname '{rulename}.{jobid}' \
	--rerun-incomplete \
	--nolock \
	--verbose \
	-k -p \
	-j 3000 \
	--stats cr_pipeline_${NOW}.stats \
	--cluster "sbatch --mail-type=FAIL -o logs/{params.rulename}.%j.o {params.batch}" \
	>& cr_pipeline_${NOW}.log

# Summary
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --summary

## DRY Run with Print out the shell commands that will be executed
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dryrun -p -r
# snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --dryrun -p -r

#DAG
 # snakemake --directory $WORK_DIR --snakefile $SNAKEFILE  --dag | dot -Tpng > dag.png

#Rulegraph
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE  -n --forceall --rulegraph | dot -Tpng > rulegraph.png

# Mail Rulegraph and DAG to self
#  echo DAG |mutt -s "DAG" -a dag.png -a rulegraph.png -- brooksma@mail.nih.gov
