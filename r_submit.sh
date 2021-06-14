#!/bin/bash

#---------------------------------------------------------------------------------
# Account information

#SBATCH --account=LoremIpsum        

#---------------------------------------------------------------------------------
# Resources requested

#SBATCH --partition=LoremIpsum     
#SBATCH --cpus-per-task=1          # number of CPUs requested (for parallel tasks)
#SBATCH --mem-per-cpu=16G          # requested memory
#SBATCH --time=2-00:00:00          # wall clock limit (d-hh:mm:ss)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=LoremIpsum@LoremIpsum.com

#---------------------------------------------------------------------------------
# Job specific name (helps organize and track progress of jobs)

#SBATCH --job-name=LoremIpsum    	        # user-defined job name

#---------------------------------------------------------------------------------
# Print some useful variables

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

#---------------------------------------------------------------------------------
# Load necessary modules for the job

module load R

#---------------------------------------------------------------------------------
# Commands to execute below...

m=$1
s=$2 
x=$3
b=$4
e=$5
g=$6
r=$7
c=$8
n=$9
shift
p=$9
shift
d=$9
shift
v=$9

echo "$m $s $x $b $e $g $r $c $n $p $d $v"
Rscript single_experiment.R -m $m -s $s -x $x -b $b -e $e -g $g -r $r -c $c -n $n -p $p -d $d -v $v