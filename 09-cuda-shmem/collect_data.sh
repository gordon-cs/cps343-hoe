#!/bin/bash

PATH=/shared/cuda-samples/bin/x86_64/linux/release:$PATH
gpus="P1000 RTX3000 RTXA5000"
sizes="08 16 32"
types="Global Shared"
dim0=1580
dim1=1620

for gpu in ${gpus}
do
    for k in ${sizes}
    do
	# clear out any previous versions of files we will build
	for t in ${types}
	do
	    rm -f $t.${gpu}.$k
	done

	# sweep through a range of matrix dimesions
	for ((n=dim0;n<=dim1;n++))
	do
	    srun -p mp,ml --gres=gpu:${gpu} ./matmul$k $n | tee tmp.$$
	    for t in ${types}
	    do
		echo "$n $(grep $t tmp.$$ | awk '{print $5}')" >> $t.${gpu}.$k
	    done
	done
    done
done

# clean up
rm -f tmp.$$
