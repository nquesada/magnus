#########################
# Makefile for Cubature #
# 24-X-13          	#
# Nicolas Quesada  	#
#########################
CC=gcc
CFLAGS=-O3 -Wall -I/usr/include -I. -I./cubature-1.0 -DHCUBATURE
LFLAGS=-lm -L/usr/lib -lgsl -lgslcblas
CUBDIR=./cubature-1.0 #Directory where the cubature library is

%.out:%.o magnusint.o c_magnusint.out magnus.o c_magnus.c functionF.o cubature-1.0/hcubature.c
	$(CC)  $^ $(LFLAGS) -o $@

c_test_fun_int.out:c_test_fun_int.o magnusint.o c_magnusint.o magnus.o c_magnus.c functionF.o cubature-1.0/hcubature.c
	$(CC)  $^ $(LFLAGS) -o $@

magnus.o: magnus.c functionF.o

c_magnus.o: c_magnus.c functionF.o

magnusint.o: magnus.o functionF.o cubature-1.0/hcubature.c

c_magnusint.o: c_magnus.o functionF.o cubature-1.0/hcubature.c


clean:
	rm -rf *~
	rm -rf *.o
	rm -rf *.out
