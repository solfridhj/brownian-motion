CXX = g++
CXXF = -g -std=c++14 -Wall
CXXFLAGS =  -O3 -larmadillo -llapack -lblas

main: main.o Particle.o 
	$(CXX) -o main main.o Particle.o  $(CXXFLAGS)

main.o: main.cpp Particle.h
	$(CXX) -c main.cpp $(CXXFLAGS)

Particle.o: Particle.h


