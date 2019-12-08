#Yasmine Brooks
#Dr. Pounds
#Numerical Methods Project 1 Makefile
#November 4, 2019

CXX = g++
CFLAGS = -std=c++11
LDLIBS = -lm -lgsl
HEADERS = globals.h constants.h prototypes.h structs.h
OBJS = integrals.o roots.o smooth.o read.o dft.o

debug ?= n
ifeq ($(debug), y)
	CFLAGS += -g -DDEBUG
else
	CFLAGS += -O2
endif

all: main tags

main : main.o $(OBJS)
	$(CXX) $(CFLAGS) -o main main.o $(OBJS) $(LDLIBS)
 
main.o : main.cpp $(HEADERS)
	$(CXX) $(CFLAGS)  main.cpp -c
        
integrals.o : integrals.cpp $(HEADERS)
	$(CXX) $(CFLAGS)  integrals.cpp -c

roots.o : roots.cpp $(HEADERS)
	$(CXX) $(CFLAGS) roots.cpp -c
  
smooth.o : smooth.cpp $(HEADERS)
	$(CXX) $(CFLAGS) smooth.cpp -c
 
read.o : read.cpp $(HEADERS)
	$(CXX) $(CFLAGS) read.cpp -c

dft.o : dft.cpp $(HEADERS)
	$(CXX) $(CFLAGS) dft.cpp -c
  
clean:
	rm *.o

pristine:
	rm *.o
	rm main
	rm tags
	touch *

tags: 
	ctags *.h *.cpp
