# File containing the main function
MAIN = t3

# Compiler
CC ?= gcc

# Name for executable
PROGRAM = $(MAIN)

# Include all your object files
OBJECTS = $(MAIN).o 

# Libraries to be included 
LIBS = \
	   -lm \
	   -lgsl \
	   -lgslcblas

# Standard run for make
CFLAGS = \
		 -O1 
		 #-Wall \ 
		 #-Werror \
		 #-pedantic
# I consider to include these flags to be good practice  

# Flags usedfor address-sanitizing debugging
CFLAGS_ADDRESS = \
	 	 -Werror \
	 	 -pedantic \
	 	 -g \
	 	 -fno-omit-frame-pointer \
	 	 -fsanitize=address

# Flags used for memory-sanitizing debugging
CFLAGS_MEMORY = \
	 	 -Werror \
		 -pedantic \
	 	 -g \
	 	 -fPIE \
	 	 -pie \
		 -fsanitize-memory-track-origins \
	 	 -fno-omit-frame-pointer \
	 	 -fsanitize=memory


# Commands for make
# Standard make
all: $(MAIN)

# Build .o files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)


# Debugging
# Adress sanitizing
address: clean
	$(CC) $(CFLAGS_ADDRESS) $(MAIN).c -o $@ $(LIBS)

# Memory sanitizing
memory: clean
	$(CC) $(CFLAGS_MEMORY) $(MAIN).c -o $@ $(LIBS)

# Light clean.
# Removes executable and objective files
clean:
	rm -rf *.o
	rm -rf $(MAIN)
	touch *.c

# Big clean
# Removes everything but header and source code
pushclean:	
	rm -rf *.o
	rm -rf $(MAIN)
	rm -rf *.csv
	rm -rf *.pdf
	touch *.c

