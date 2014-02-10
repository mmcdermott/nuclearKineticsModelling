CXX = g++
WARNINGFLAGS = -pedantic -Wall
CXXFLAGS = -std=c++11 -O2 $(WARNINGFLAGS)

TARGETS = main

all: $(TARGETS)

main: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean: 
	rm -f *.o main
