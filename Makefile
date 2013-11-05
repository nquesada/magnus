#########################
# Makefile for Cubature #
# 24-X-13          	#
# Nicolas Quesada  	#
#########################
CC=gcc
CFLAGS=-O3 -Wall -I/usr/include -I. -I./cubature-1.0 -DHCUBATURE
LFLAGS=-lm -L/usr/lib -lgsl -lgslcblas
CUBDIR=./cubature-1.0 #Directory where the cubature library is

%.out:%.o magnus.o functionF.o cubature-1.0/hcubature.c
	$(CC)  $^ $(LFLAGS) -o $@


magnus.o: magnus.c functionF.o



#test.out:test.o magnus.o functionF.o $(CUBDIR)/hcubature.c
#	$(CC)  $^ $(LFLAGS) -o $@

#nico.out: nico.o

clean:
	rm -rf *~
	rm -rf *.o
	rm -rf *.out