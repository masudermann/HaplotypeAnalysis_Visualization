# A script for parallelizing the clustering analysis

log="A"
window=100000
step=25000
cutoff=10
dmin=2
dmax=100
dstep=4
nohup python3 -u cluster_haplotypes.py SL4.0ch00.vcf ch00 $window $step $cutoff $dmin $dmax $dstep > ch00${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch01.vcf ch01 $window $step $cutoff $dmin $dmax $dstep > ch01${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch02.vcf ch02 $window $step $cutoff $dmin $dmax $dstep > ch02${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch03.vcf ch03 $window $step $cutoff $dmin $dmax $dstep > ch03${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch04.vcf ch04 $window $step $cutoff $dmin $dmax $dstep > ch04${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch05.vcf ch05 $window $step $cutoff $dmin $dmax $dstep > ch05${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch06.vcf ch06 $window $step $cutoff $dmin $dmax $dstep > ch06${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch07.vcf ch07 $window $step $cutoff $dmin $dmax $dstep > ch07${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch08.vcf ch08 $window $step $cutoff $dmin $dmax $dstep > ch08${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch09.vcf ch09 $window $step $cutoff $dmin $dmax $dstep > ch09${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch10.vcf ch10 $window $step $cutoff $dmin $dmax $dstep > ch10${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch11.vcf ch11 $window $step $cutoff $dmin $dmax $dstep > ch11${log}.log &
nohup python3 -u cluster_haplotypes.py SL4.0ch12.vcf ch12 $window $step $cutoff $dmin $dmax $dstep > ch12${log}.log &