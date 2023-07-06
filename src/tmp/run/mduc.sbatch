#!/bin/bash

## Job name and files
#SBATCH --job-name=mduc

## OUTPUT AND ERROR
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err

## Initial working directory
##SBATCH -D $PWD

## specify your mail address (send feedback when job is done)
##SBATCH --mail-type=begin  				 # send email when process begins...
##SBATCH --mail-type=end    				 # ...and when it ends...
##SBATCH --mail-type=fail   		 	         # ...or when it fails.
##SBATCH --mail-user=luca.toni@itv.rwth-aachen.de	 # send notifications to this emai REPLACE WITH YOUR EMAIL

## Setup of execution environment 
#SBATCH --export=NONE

## account and partition selection
#SBATCH --account='itv'
#SBATCH --partition='ih'

## specific partition selection
#SBATCH -C 'broadwell'
##SBATCH -w 'lnih021 lnih022 lnih023 lnih024 lnih025'
##SBATCH -w 'lnih026 lnih027 lnih028 lnih029 lnih030'

###  for citv 1-4 ###
##SBATCH -C 'westmere'

###  for citv 5   ###
##SBATCH -C 'ivybridge'

## Request the number of nodes
#SBATCH --nodes=5           # number of nodes

##### ITV clusters have 12 cores
## Set tasks per node
#SBATCH --tasks-per-node=24    # NOTE new machine has 48 core per node but for our partition you should set 12 

## Set cpus per task (multiple tasks for hyperthreading)
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1


## memory per cpu
##SBATCH --mem-per-cpu=7600


## do not share nodes 
##SBATCH --exclusive

## OPTIONAL
## set max. file size to 500Gbyte per file
##SBATCH -F 500000000

## execution time in [hours]:[minutes]:[seconds]
## recommended: less than 24:00
#SBATCH --time=1:00:00

############################################### 
###############################################
########## END OF SLURM INSTRUCTION ###########
##!!! DO NOT PLACE SHELL VARIABLES ABOVE!!!####
###############################################
###############################################

## Optional lines to remove all modules loaded and load the one disired
## Not need to use them if you use the default modules
## In the case the code has been compiled with the intel2016, change the following line accordingly 

#module purge; module load DEVELOP DEPRECATED
#module load intel
#module load intelmpi/2018.3.222


#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/itv/compile_centos7/lib/sundials/2.7.0_intel19.0_intelmpi_2018_mkl/Install_DIR/lib/
#

## #executable name (e.g. arts_cf)
exe='mduc'

### executable arguments
args='-i mduc.in'

date
mpiexec.hydra $FLAGS_MPI_BATCH ./$exe $args

date

