include Makefile.inc
strain: main.cc *h
	$(CC) -o strain main.cc $(LDFLAGS) $(DEBUGFLAGS)

test: test.cc *h
	g++ -std=gnu++0x -o test test.cc -O2
	./test

run:	strain
	./strain $(SEED)

ps:	
	head -1 tree > _sorted
	gawk "{if(NR!=1) print $0}" tree | sort -r >> _sorted
	./tree2ps _sorted > test.ps
	rm _sorted

	ps2eps -f test.ps
	gv test.eps

bak:
	mv -f single0* out logtime tree trash/

	
#1321731977
