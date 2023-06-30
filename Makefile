CC=gcc -m64 
CFLAGS=-c -Wall -lm -ldl

all: Mu2P

Mu2P: main.o matrix_element.o brody.o potential.o angular.o romberg.o
	$(CC) main.o matrix_element.o brody.o potential.o angular.o romberg.o -o Mu2P -lm -ldl -lgsl -lgslcblas

main.o: main.c
	$(CC) $(CFLAGS) main.c

matrix_element.o: matrix_element.c
	$(CC) $(CFLAGS) matrix_element.c

potential.o: potential.c
	$(CC) $(CFLAGS) potential.c

angular.o: angular.c
	$(CC) $(CFLAGS) angular.c

brody.o: brody.c
	$(CC) $(CFLAGS) brody.c

romberg.o: romberg.c
	$(CC) $(CFLAGS) romberg.c

clean:
	rm -rf *.o Mu2P
