#!/bin/bash

# No permutations for heteroskedastic errors

# m: method
for m in 'blpr' 'dlasso' 'silm'; do
# for m in 'silm'; do
# s: sparsity
for s in 4 15; do
# x: x_design
for x in 'N1' 'G1' 'N2' 'TG' 'TGM' 'WB'; do
# b: b_design
for b in 'D1'; do
# e: e_design
for e in 'N1' 'G1' 'N2' 'WB'; do
# g: perm_design
for g in 'perm'; do
# r: # sim
for r in 200; do
# c: seed index
for c in 0 200 400 600 800; do
# for c in 0; do
# n: # obs
for n in 100; do
# p: # dim
for p in 300; do
# d: # draws
for d in 1000; do
# v: standardize
for v in 1; do
    # sbatch --account=pi-mkolar rcc_r_submit.sh $m $s $x $b $e $g $r $c $n $p $d $v &
    sbatch --account=pi-mkolar r_submit.sh $m $s $x $b $e $g $r $c $n $p $d $v &
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