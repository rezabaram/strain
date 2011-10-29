all:	strain
	./strain > out 
strain: main.cc *h
	g++ -o strain main.cc -O2
