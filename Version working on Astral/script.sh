#!/bin/bash
##
## MPI submission script for PBS ASTRAL
## ------------------------------------
##
## Follow the 4 steps below to configure. If you edit this from Windows,
## *before* submitting via "qsub" run "dos2unix" on this file - or you will
## get strange errors. You have been warned.
##
## STEP 1:
## The following line contains the job name:
##
#PBS -N 16procs
##
## STEP 2:
##
## The select line below selects 1 chunk of 16 cpus
## Maximum values for ppn is 16.
##
#PBS -l nodes=1:ppn=16
##Course Announcement →  MSc in Computational & Software Techniques in Engineering 2016/2017Internships and Thesis Open Dismiss


## STEP 3:
##
## Select correct queue:
##
## express    - 2 hours
## short      - 8 hours
## medium     - 1 day
## long       - 3 days
## very_long  - 5 days
## extra_long - 10 days (by special arrangement)
## cismg      - 32 cpu nodes (by special arrangement)
##
#PBS -q express
##
## STEP 4:
##
##
## DO NOT CHANGE the following lines
##------------------------------------------------
#PBS -k oe
##
## Change to working directory
outfl=/panfs/storage/s261745/Assignment/a.txt
cd $PBS_O_WORKDIR
##
## Find hosts
cat $PBS_NODEFILE > hosts
##
## Set up INTEL environment.
. /etc/profile.d/modules.sh
module load icc
module load ifort
module load impi
. iccvars.sh intel64
. ifortvars.sh intel64
. mpivars.sh
##
##-------------------------------------------------
##
## STEP 5:
##
## Put correct parameters in mpirun execution line
## below:
##
mpirun -f hosts -np 16 ./MPI_CMP 2>&1 $outfl
