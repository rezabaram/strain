strain: main.cc *h
	g++-fsf-4.5 -std=gnu++0x -g -o strain main.cc -O2

test: test.cc *h
	g++-fsf-4.5 -std=gnu++0x -o test test.cc -O2
	./test

run:	strain
	./strain 
