CC=gcc -pthread -lm
CFLAGS=-I. -O2
DEPS=

all: newton

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

newton: main.o
	$(CC) -o newton main.o

clean:
	rm -f *.o $(objects) newton
