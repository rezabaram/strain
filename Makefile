include Makefile.inc
strain: main.cc *h
	$(CC) -o strain main.cc $(LDFLAGS) $(DEBUGFLAGS)

test: test.cc *h
	g++ -std=gnu++0x -o test test.cc -O2
	./test

run:	strain
	./strain $(SEED)

ps:	
	head -1 $(FILE) > _sorted
	gawk "{if(NR!=1) print $0}" $(FILE) | sort -r >> _sorted
	./genps _sorted > $(FILE).ps
	rm _sorted
	ps2eps -f $(FILE).ps
	gv $(FILE).eps

bak:
	rm -rf trash/*
	mv -f tree0* single0* out logtime tree trash/

clean:
	rm -f tree0* single0* out logtime tree 

	
#1321731977
