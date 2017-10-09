ifeq ($(ARCH),LCS_KNC) 
	CFLAGS=-g -c -openmp  -O3 -DLCS_KNC
	LDFLAGS=-openmp
	ICC=mpiicpc
	SOURCES = main.c node_compute.c global.c lcs.c seqs_handle.c file.c 
	OBJS= node_compute.o main.o global.o lcs.o seqs_handle.o file.o 
else
	CFLAGS=-g -c -qopenmp  -O3 -DLCS_KNL
	LDFLAGS=-qopenmp 
	ICC=mpiicpc
	SOURCES = main.c node_compute.c global.c seqs_handle.c file.c lcs_knl.c
	OBJS= node_compute.o main.o global.o seqs_handle.o file.o lcs_knl.o
endif


main: $(OBJS)
	$(ICC) $(OBJS) $(LDFLAGS) -o main

node_compute.o: node_compute.c node_compute.h
	$(ICC) $(CFLAGS) node_compute.c -o node_compute.o

global.o:global.c global.h
	$(ICC) $(CFLAGS) global.c -o global.o

lcs.o:lcs.c lcs.h
	$(ICC) $(CFLAGS) lcs.c -o lcs.o

ifeq ($(ARCH),LCS_KNL)
lcs_knl.o:lcs_knl.c lcs_knl.h
	$(ICC) $(CFLAGS) lcs_knl.c
endif

seqs_handle.o:seqs_handle.c seqs_handle.h
	$(ICC) $(CFLAGS) seqs_handle.c -o seqs_handle.o 

file.o:file.c file.h
	$(ICC) $(CFLAGS) file.c -o file.o

main.o: main.c
	$(ICC) $(CFLAGS) main.c -o main.o

clean:
	rm *.o main
