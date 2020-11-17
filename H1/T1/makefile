
CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm

HEADERS = H1lattice.h H1potential.h write_to_file.h
OBJECTS = H1lattice.o H1potential.o write_to_file.c H1main.o
PROGRAM = MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	rm $(PROGRAM)
	touch *.c

pushclean:
	rm -f *.o
	rm $(PROGRAM)
	rm -rf *.csv
	rm -rf *.pdf
	touch *.c


