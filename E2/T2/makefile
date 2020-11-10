main = t2

CC ?= gcc # c-compiler
CFLAGS = # compiler  flags  such as -O3
LIBS = -lm -lgsl -lgslcblas  # libraries  such as -lm -lgsl -lgslcblas
OBJECTS = $(main).o fft.o # object  files  to  produce
PROGRAM = $(main) # name of  executable  to be  built# build .o files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
	# compile  program
all: $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -rf *.o
	rm -rf $(main)
	touch *.c

pushclean:	
	rm -rf *.o
	rm -rf $(main)
	rm -rf *.csv
	rm -rf *.pdf
	touch *.c
	

