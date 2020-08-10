#!/bin/bash

# No sign flips for gamma errors

# s: sparsity
for s in 3 10; do
# x: x_design
# for x in 'N1' 'G1' 'N2' 'TG' 'TGM'; do
for x in 'WB' 'TWB'; do
# b: b_design
for b in 'D1'; do
# e: e_design
for e in 'N1' 'N2' 'HG' 'HMG'; do
# g: perm_design
for g in 'sign'; do
# r: # sim
for r in 50; do
# n: # obs
for n in 50; do
# p: # dim
for p in 100; do
# d: # draws
for d in 1000; do
# i: # solves
for i in 500; do
# m: rescale errors
for m in 1; do
  sbatch --account=pi-mkolar r_submit.sh $s $x $b $e $g $r $n $p $d $i $m &
done;
done;
done;
done;
done;
done;
done;
done;
done;
done;
done;