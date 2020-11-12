# File containing the main function
MAIN = t2

# Compiler
CC ?= gcc

# Name for executable
PROGRAM = $(MAIN)

# Include all your object files
OBJECTS = $(MAIN).o fft.o 

# Libraries to be included 
LIBS = -lm \ 		# Math		
	   -lgsl \		# GSL
	   -lgslcblas   # BLAS for GSL

# Standard run for make
CFLAGS = -O1 \		# Optimizing
		 -Wall \	# Enables all errors
		 -Werror \	# Makes warnings into error
		 -pedantic	# Code needs to comply with ISO C
# I consider to include these flags to be good practice  

# Flags usedfor address-sanitizing debugging
CFLAGS_ADDRESS = \
	 	 -Werror \					# Makes warning into errors
	 	 -pedantic \				# Code need to comply with ISO C
	 	 -g \						# Makes stack trace human readable
	 	 -fno-omit-frame-pointer \	# Nicer printed stack trace
	 	 -fsanitize=address			# Enables address sanitizer

# Flags used for memory-sanitizing debugging
CFLAGS_MEMORY = \								
	 	 -Werror \							# Makes warning into erros
	 	 -pedantic \						# Code needs to comply with ISO C
	 	 -g \								# Makes stack trace human readable
	 	 -fPIE \							# Needed for -pie
	 	 -pie \								# Generates position depedent code
		 -fsanitize-memory-track-origins \	# Memory sanitizer tracks whole path of uninitalized value
	 	 -fno-omit-frame-pointer \			# Nicer printed stack trace
	 	 -fsanitize=memory					# Enables memory sanitizer


# Definition of make commands
# Standard make
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
	
all: $(MAIN)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

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
	
# Debugging
# Adress sanitizing
address: clean
	$(CC) $(CFLAGS_ADDRESS) $(MAIN).c -o $@ $(LIBS)

# Memory sanitizing
memory: clean
	$(CC) $(CFLAGS_MEMORY) $(MAIN).c -o $@ $(LIBS)

