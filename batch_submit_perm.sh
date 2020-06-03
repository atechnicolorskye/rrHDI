#!/bin/bash

# d: sparsity
for d in 3 10; do
# for d in 3; do
# x: x_design
for x in 'N1' 'G1' 'N2' 'TG' 'TGM'; do
# for x in 'N1'; do
# b: b_design
for b in 'D1'; do
# e: e_design
for e in 'N1' 'G1' 'N2' 'HG' 'HMG'; do
# for e in 'N1'; do
# g: perm_design
for g in 'perm'; do
# r: # sim
for r in 50; do
# for r in 10; do
# # n: # obs
# for n in 50; do
# # p: # dim
# for p in 200 300 400; do
# s: # solves
for s in 500; do
  # echo "$d $x $b $e $g $r $n $p $s"
  sbatch --account=pi-mkolar r_submit.sh $d $x $b $e $g $r 0 0 $s &
# done;
# done;
done;
done;
done;
done;
done;
done;
done;

# no sign flips for gamma