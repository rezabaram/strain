strain: main.cc *h
	g++  -W -std=gnu++0x -g -o strain main.cc -O2

test: test.cc *h
	g++ -std=gnu++0x -o test test.cc -O2
	./test

run:	strain
	./strain $(SEED)

ps:	tree
	./tree2ps tree > test.ps
	ps2eps -f test.ps
	gv test.eps

bak:
	mv -f single0* out logtime tree trash/

	
#1321731977
