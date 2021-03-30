#!/bin/bash

# No permutations for heteroskedastic errors

# s: sparsity
for s in 3 10; do
# x: x_design
for x in 'WB'; do
# b: b_design
for b in 'D1'; do
# e: e_design
for e in 'N1' 'G1' 'N2' 'WB'; do
# g: perm_design
for g in 'perm'; do
# r: # sim
for r in 100; do
# c: seed index
for c in 0 100 200 300 400 500 600 700 800 900; do
# n: # obs
for n in 50; do
# p: # dim
for p in 100; do
# d: # draws
for d in 1000; do
    sbatch --account=pi-mkolar m_rr_hdi_r_submit.sh $s $x $b $e $g $r $c $n $p $d &
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

# s: sparsity
for s in 3 10; do
# x: x_design
for x in 'N1' 'G1' 'N2' 'TG' 'TGM'; do
# b: b_design
for b in 'D1'; do
# e: e_design
for e in 'WB'; do
# g: perm_design
for g in 'perm'; do
# r: # sim
for r in 100; do
# c: seed index
for c in 0 100 200 300 400 500 600 700 800 900; do
# n: # obs
for n in 50; do
# p: # dim
for p in 100; do
# d: # draws
for d in 1000; do
    sbatch --account=pi-mkolar m_rr_hdi_r_submit.sh $s $x $b $e $g $r $c $n $p $d &
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