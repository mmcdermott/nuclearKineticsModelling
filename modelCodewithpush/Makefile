CXX = g++
WARNINGFLAGS = -pedantic -Wall
CXXFLAGS = -std=c++11 -O2 $(WARNINGFLAGS)

TARGETS = clean mtKineticModel

all: $(TARGETS)

mtKineticModel: Vector.o main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

vectorTest: Vector.o vectorTest.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean: 
	rm -f *.o mtKineticModel
