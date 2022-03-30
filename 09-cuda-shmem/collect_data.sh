#!/bin/bash

PATH=/shared/cuda-samples/bin/x86_64/linux/release:$PATH
arch=$(deviceQuery | grep "Device 0" | cut -d: -f2 | sed -e 's/[ "]//g')
types="Global Shared"
sizes="08 16 32"
dim0=1580
dim1=1620

for k in ${sizes}
do
    # clear out any previous versions of files we will build
    for t in ${types}
    do
	rm -f $t.${arch}.$k
    done

    # sweep through a range of matrix dimesions
    for ((n=dim0;n<=dim1;n++))
    do
	./matmul$k $n | tee tmp.$$
	for t in ${types}
	do
	    echo "$n $(grep $t tmp.$$ | awk '{print $5}')" >> $t.${arch}.$k
	done
    done
done

# clean up
rm -f tmp.$$
