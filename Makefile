all:
	cd src && make
	cp src/matchclips .
	cp src/cnvtable .

clean:
	cd src && make clean
	rm matchclips cnvtable
