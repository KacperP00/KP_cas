#!/bin/zsh

## Job name and files
#SBATCH --job-name=test_job

## OUTPUT AND ERROR
#SBATCH -o job.%j.out
#SBATCH -o job.%j.out

## Initial working directory
#SBATCH -D .

## specify your mail address (send feedback when job is done)
#SBATCH --mail-type=begin  				 # send email when process begins...
#SBATCH --mail-type=end    				 # ...and when it ends...
#SBATCH --mail-type=fail   		 	         # ...or when it fails.
##SBATCH --mail-user=<YOUREMAIL>	 # send notifications to this emai REPLACE WITH YOUR EMAIL

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


## Set tasks per node
#SBATCH --tasks-per-node=1    # NOTE new machine has 48 core per node


## Set cpus per task (multiple tasks for hyperthreading)
#SBATCH --cpus-per-task=1


## memory per cpu
##SBATCH --mem-per-cpu=7600


## do not share nodes 
##SBATCH --exclusive

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

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/itv/compile_centos7/lib/sundials/2.7.0_intel19.0_intelmpi_2018_mkl/Install_DIR/lib/

## #executable name (e.g. arts_cf)
exe='cas_opt'

### executable arguments
args='cas.in'

date
$exe $args
date
