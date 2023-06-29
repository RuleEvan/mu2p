CC=gcc -m64 
CFLAGS=-c -Wall -lm -ldl

all: Mu2E

Mu2E: main.o matrix_element.o wfn.o charge.o file_io.o multipole.o brody.o potential.o angular.o romberg.o
	$(CC) main.o matrix_element.o wfn.o charge.o file_io.o multipole.o brody.o potential.o angular.o romberg.o -o Mu2E -lm -ldl -lgsl -lgslcblas

main.o: main.c
	$(CC) $(CFLAGS) main.c

matrix_element.o: matrix_element.c
	$(CC) $(CFLAGS) matrix_element.c

wfn.o: wfn.c
	$(CC) $(CFLAGS) wfn.c

harmonic.o: harmonic.c
	$(CC) $(CFLAGS) harmonic.c

potential.o: potential.c
	$(CC) $(CFLAGS) potential.c

charge.o: charge.c
	$(CC) $(CFLAGS) charge.c

file_io.o: file_io.c
	$(CC) $(CFLAGS) file_io.c

multipole.o: multipole.c
	$(CC) $(CFLAGS) multipole.c

angular.o: angular.c
	$(CC) $(CFLAGS) angular.c

brody.o: brody.c
	$(CC) $(CFLAGS) brody.c

romberg.o: romberg.c
	$(CC) $(CFLAGS) romberg.c

clean:
	rm -rf *.o Mu2E
