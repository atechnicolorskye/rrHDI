#!/bin/bash

#---------------------------------------------------------------------------------
# Account information

#SBATCH --account=staff              # basic (default), staff, phd, faculty

#---------------------------------------------------------------------------------
# Resources requested

#SBATCH --partition=standard       # standard (default), long, gpu, mpi, highmem
#SBATCH --cpus-per-task=1          # number of CPUs requested (for parallel tasks)
#SBATCH --mem-per-cpu=16G          # requested memory
#SBATCH --time=2-00:00:00          # wall clock limit (d-hh:mm:ss)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sikai.lee@chicagobooth.edu

#---------------------------------------------------------------------------------
# Job specific name (helps organize and track progress of jobs)

#SBATCH --job-name=m_rr_hdi    	  # user-defined job name

#---------------------------------------------------------------------------------
# Print some useful variables

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

#---------------------------------------------------------------------------------
# Load necessary modules for the job

module purge
module load R/3.6/3.6.2
module load gcc/9.2.0

#---------------------------------------------------------------------------------
# Commands to execute below...

s=$1 
x=$2
b=$3
e=$4
g=$5
r=$6
c=$7
n=$8
p=$9
shift
d=$9

echo "$s $x $b $e $g $r $c $n $p $d"
Rscript single_experiment_m_rr_hdi.R -s $s -x $x -b $b -e $e -g $g -r $r -c $c -n $n -p $p -d $d
