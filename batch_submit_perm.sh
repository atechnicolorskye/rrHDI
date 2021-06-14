#!/bin/bash

# No permutations for heteroskedastic errors

# m: method
for m in 'blpr' 'dlasso' 'hdi' 'silm' 'rr'; do
# s: sparsity
for s in 4 15; do
# x: x_design
for x in 'N1' 'G1' 'N2' 'NT' 'GT' 'WB'; do
# b: b_design
for b in 'D1'; do
# e: e_design
for e in 'N1' 'G1' 'N2' 'WB'; do
# g: perm_design
for g in 'perm'; do
# r: # sim
for r in 100; do
# # n: # obs
for n in 50; do
# p: # dim
for p in 100; do
# d: # draws
for d in 1000; do
# v: standardize
for v in 0 1; do
    sbatch r_submit.sh $m $s $x $b $e $g $r $n $p $d $v &
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
done;