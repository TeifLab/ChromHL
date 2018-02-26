#!/bin/sh

# SGE job script for sumbitting matlab jobs to an SGE cluster queue.
# This is submitted to the queue by the job.sh script. It simply runs
#  matlab with the correct arguments.
# By David Black-Schaffer, June 2007.
# Permission to use and modify this script is granted.
# I am not responsible for any errors in this script, so be forewarned!


#$ -j n
#$ -o job-nobackup.$JOB_ID.$TASK_ID.out
#$ -e job-nobackup.$JOB_ID.$TASK_ID.err
#$ -cwd
#$ -m beas
#$ -M gthorn@essex.ac.uk
#$ -t 1:37
#$ -tc 20

echo "Starting job: $SGE_TASK_ID"

# Modify this to use the path to matlab for your system

matlab -nojvm -nodisplay -r multiple_run_mcd_20  < /dev/null

echo "Done with job: $SGE_TASK_ID"

