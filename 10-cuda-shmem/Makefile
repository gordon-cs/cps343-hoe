matmul:	matmul.cu wtime.o wtime.h
	nvcc -arch=sm_30 -O2 -o matmul matmul.cu wtime.o -lrt

wtime.o: wtime.c wtime.h
	gcc -c -O2 -o wtime.o wtime.c

clean clobber:
	rm -f matmul wtime.o
