#!/bin/bash

arch=$(deviceQuery | grep "Device 0" | cut -d: -f2 | sed -e 's/[ "]//g')

for k in 08 16 32
do
    for ((n=1580;n<=1620;n++))
    do
	./matmul$k $n
    done | tee data.${arch}.$k
    for t in Global Shared
    do
	grep $t data.${arch}.$k | awk '{print $5}' > $t.${arch}.$k
    done
done
