CC=gcc
CFLAGS=-c -Wall
MATH=-lm
OBJECTS=algorithm.o mersenne.o sosfunctions.o

all: alg

alg: $(OBJECTS)
	$(CC) $(OBJECTS) -o alg $(MATH)

algorithm.o: algorithm.c
	$(CC) $(CFLAGS) algorithm.c

mersenne.o: mersenne.c
	$(CC) $(CFLAGS) mersenne.c $(MATH)

sosfunctions.o: sosfunctions.c
	$(CC) $(CFLAGS) sosfunctions.c $(MATH)

clean: 
	rm $(OBJECTS)

