all:
	cd src && make
	cp src/matchclips .

clean:
	cd src && make clean
	rm matchclips
