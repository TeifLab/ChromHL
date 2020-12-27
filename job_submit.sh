#!/bin/sh

# This script is for submitting multiple runs to the an SGE cluster queue
# in MATLAB.
#
# Alter the -m for notification options
#           -M for email notification
#           -t for number of subtasks (passed as env variable SGE_TASK_ID)
#           -tc for number of concurrent tasks to run
#
# 
# Original header below:
#
## SGE job script for sumbitting matlab jobs to an SGE cluster queue.
## This is submitted to the queue by the job.sh script. It simply runs
##  matlab with the correct arguments.
## By David Black-Schaffer, June 2007.
## Permission to use and modify this script is granted.
## I am not responsible for any errors in this script, so be forewarned!


#$ -j n
#$ -o job-nobackup.$JOB_ID.$TASK_ID.out
#$ -e job-nobackup.$JOB_ID.$TASK_ID.err
#$ -cwd
#$ -m beas
#$ -M <email> for notifications
#$ -t 1:37
#$ -tc 20

echo "Starting job: $SGE_TASK_ID"

# Run matlab without display, without jvm 
# -r is the option for which script to run, in this case
# multiple_run_chromhl (do not add .m)
# "< /dev/null" prevents MATLAB giving errors running
# non-interactively

matlab -nojvm -nodisplay -r multiple_run_chromhl  < /dev/null

echo "Done with job: $SGE_TASK_ID"

# Note that the .out file will have
#
# Warning: no access to tty (Bad file descriptor).
# Thus no job control in this shell.
#
# at top - this can be ignored safely.
