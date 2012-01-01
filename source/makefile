#Type "make" to compile all then "./multirun"
#Type "make simulation" to compile only the main simulation program

all : multirun

simulation: main.o space.o box.o
	g++ main.o space.o box.o -o simulation
	
box.o: box.cpp box.h
	g++ -c box.cpp
	
main.o: main.cpp space.h box.h
	g++ -c main.cpp
	
space.o: space.cpp space.h
	g++ -c space.cpp
	
multirun: multirun.cpp simulation
	g++ multirun.cpp -o multirun
