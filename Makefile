#########################
# Makefile for Cubature #
# 24-X-13          	#
# Nicolas Quesada  	#
#########################
CC=gcc
CFLAGS=-O3 -Wall -I/usr/include -I. -I./cubature-1.0
LFLAGS=-lm -L/usr/lib -lgsl -lgslcblas
CUBDIR=./cubature-1.0 #Directory where the cubature library is

magnus.o: magnus.c functionF.o

%.out:%.o magnus.o functionF.o $(CUBDIR)/hcubature.c
	$(CC)  $^ $(LFLAGS) -o $@

#nico.out: nico.o

clean:
	rm -rf *~
	rm -rf *.o
	rm -rf *.out