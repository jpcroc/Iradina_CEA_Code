# This makefile can be used to compile iradina.

# For Linux or UNIX systems it is recommend to use the gcc for compilation
# On Windows, using mingw32 is recommended
#Default:
#CC = x86_64-w64-mingw32-gcc
#select compiler by OS:
#OS := $(shell uname -s)
#ifeq ($(OS),Linux)
  CC=gcc
#else
  #assume windows:
#  CC=x86_64-w64-mingw32-gcc
#endif

# Add -m32 for 32 bit compilation on 64 bit systems:
CFLAGS = -O1 -Wall -fcommon
LDFLAGS = -lm

########################################
# Notes on warning level:
########################################
# Using -Wall and -pedantic will return many warnings, because of non-allowed comment styles in the codes.
# I recommend using just -Wall. 

########################################
# Notes on the compiler optimization options -OX:
########################################
# It is not recommended to use -O2 or -O3. On Linux systems no problems using these options have
# been observed so far; however, on windows systems porgram crashes and hang-ups did occur when
# using -O2 or -O3. The problems occurred in modules "transport.c" and "fromcorteo.c". 
# In windows, the inverse sqr tables will not work with -O2.
# If unsure, compile all with -O1 only.

########################################
# Notes on 32/64 bit systems:
# Since iradina version 1.0.6 native 64 bit compilation is possible.
# If you want to build a 32bit version using a 64bit compiler, active the flag -m32 for the compiler and the linker.
########################################

iradina : iradina.o target.o fileio.o utils.o transport.o geometry.o fromcorteo.o
	$(CC) -o iradina iradina.o target.o fileio.o utils.o transport.o geometry.o fromcorteo.o $(LDFLAGS)

utils.o: utils.h utils.c

target.o: target.h target.c

iradina.o: iradina.h iradina.c

fileio.o: fileio.h fileio.c

transport.o: transport.h transport.c

geometry.o: geometry.h geometry.c

fromcorteo.o: fromcorteo.h fromcorteo.c
#	gcc -lm -c -O1 -Wall fromcorteo.c fromcorteo.h

.c.o :
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f iradina *.o
