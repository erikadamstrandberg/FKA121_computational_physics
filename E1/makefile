CC ?= gcc # c-compiler
CFLAGS = # compiler  flags  such as -O3
LIBS = -lm -lgsl -lgslcblas  # libraries  such as -lm -lgsl -lgslcblas
OBJECTS = E1code3.o fft.o # object  files  to  produce
PROGRAM = run # name of  executable  to be  built# build .o files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
	# compile  program
all: $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf *.o
	rm -rf *.exe
	touch *.c
