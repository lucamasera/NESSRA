CC = gcc
CXX = g++
CFLAGS = -c -O3 -Wall -Wextra -pedantic
LDFLAGS = 
CXXSOURCES = Graph.cpp utility.cpp cpdag.cpp skeleton.cpp pc.cpp
CSOURCES = erf.c
OBJECTS = $(CXXSOURCES:.cpp=.o) $(CSOURCES:.c=.o)
EXECUTABLE = pc++

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o bin/$@

.cpp.o:
	$(CXX) $(CFLAGS) $< -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf bin/$(EXECUTABLE) *.o *~

.PHONY: all clean
