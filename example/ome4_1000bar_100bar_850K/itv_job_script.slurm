#!/usr/local_rwth/bin/zsh

## Job name and files
#SBATCH --job-name=CATS_OME1

## OUTPUT AND ERROR
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err

## Initial working directory
#SBATCH -D .

## specify your mail address (send feedback when job is done)
##SBATCH --mail-type=begin  				 # send email when process begins...
#SBATCH --mail-type=end    				 # ...and when it ends...
##SBATCH --mail-type=fail   		 	         # ...or when it fails.
#SBATCH --mail-user=a.deshmukh@itv.rwth-aachen.de	 # send notifications to this emai REPLACE WITH YOUR EMAIL

## Setup of execution environment 
#SBATCH --export=NONE

## choose account
#SBATCH --account='itv'
#SBATCH --partition='ih'
##SBATCH -C 'westmere'  ### CITV1-4
##SBATCH -C 'ivybridge' ### CITV5
#SBATCH -C 'c2'	       ### CITV6-7

## Request the number of nodes
#SBATCH --nodes=1           # number of nodes

####SBATCH --ntasks=96

## Set tasks per node
#SBATCH --tasks-per-node=1    # NOTE new machine has 48 core per node


## Set cpus per task (multiple tasks for hyperthreading)
#SBATCH --cpus-per-task=1


## memory per cpu
#SBATCH --mem-per-cpu=3600


## do not share nodes 
#SBATCH --exclusive

## OPTIONAL
## set max. file size to 500Gbyte per file
##SBATCH -F 500000000

## execution time in [hours]:[minutes]:[seconds]
## recommended: less than 24:00
#SBATCH --time=24:00:00

############################################### 
###############################################
########## END OF SLURM INSTRUCTION ###########
##!!! DO NOT PLACE SHELL VARIABLES ABOVE!!!####
###############################################
###############################################

#export I_MPI_PIN=1
#export I_MPI_PIN_DOMAIN=core

## #executable name (e.g. arts_cf)
exe='cas_opt'
#exe='cas'

### executable arguments
args='cas.in'

# run executable
date
#$MPIEXEC $FLAGS_MPI_BATCH $exe $args 2>&1 | tee out.job_$SLURM_JOB_ID
$exe $args 2>&1 | tee out.job_$SLURM_JOB_ID
date
