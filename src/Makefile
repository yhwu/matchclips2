#PROGRAMS = matchclips simbp tagcnv cnvview
PROGRAMS = matchclips
SAMTOOLS = samtools-0.1.19
all: $(PROGRAMS)
TITLE = matchclip2
CC = g++
CFLAGS = -g -O2  -Wall -L/usr/lib64
#CFLAGS = -g -Wall -L/usr/lib64
INC=-I${SAMTOOLS}
LIBS =  -L${SAMTOOLS} -lbam -pthread -lm -lz
BACKUPDIR = ~/projects/$(strip $(TITLE))
BACKUPLOG = $(strip $(BACKUPDIR)/$(TITLE)).log
STAMP = $(strip $(shell  date +'%Y.%m.%d-%H.%M.%S'))
BACKUPFOLDER = $(BACKUPDIR)/$(STAMP)	

TAGCXX =  tagcnvmain.cpp tagcnv.cpp readref.cpp read_CIGAR.cpp functions.cpp
TAGHDR = $(TAGCXX:.cpp=.h)	
TAGOBJ = $(TAGCXX:.cpp=.o)	
tagcnv : $(TAGOBJ) $(TAGCXX) $(TAGHDR) Makefile
	$(CC) $(CFLAGS) $(TAGOBJ)  -o $@

MATCHCXX =  matchreadsmain.cpp matchreads.cpp preprocess.cpp exhaustive.cpp pairguide.cpp  samfunctions.cpp readref.cpp functions.cpp
MATCHHDR = color.h $(MATCHCXX:.cpp=.h)	
MATCHOBJ = $(MATCHCXX:.cpp=.o)	
matchclips : $(MATCHOBJ) $(MATCHCXX) $(MATCHHDR) Makefile ./${SAMTOOLS}/libbam.a
	$(CC) $(CFLAGS) $(MATCHOBJ) $(INC) $(LIBS) -o $@

time:
	date -u "+%a %b %d %H:%M:%S %Y" > UPDATED

CNVVIEWCXX =  cnvview.cpp filters.cpp samfunctions.cpp readref.cpp functions.cpp
CNVVIEWHDR = matchreads.h $(CNVVIEWCXX:.cpp=.h)	
CNVVIEWOBJ = $(CNVVIEWCXX:.cpp=.o)	
cnvview : $(CNVVIEWOBJ) $(CNVVIEWCXX) $(CNVVIEWHDR) Makefile ./${SAMTOOLS}/libbam.a
	$(CC) $(CFLAGS) $(CNVVIEWOBJ) $(INC) $(LIBS) -o $@

./${SAMTOOLS}/libbam.a :
	cd ${SAMTOOLS} && make libbam.a

SIMBPCXX =  simbp.cpp readref.cpp 
SIMBPHDR = $(SIMBPCXX:.cpp=.h)	
SIMBPOBJ = $(SIMBPCXX:.cpp=.o)	wgsim.o
simbp : $(SIMBPOBJ) $(SIMBPCXX) $(SIMBPHDR) wgsim.c kseq.h Makefile
	$(CC) $(CFLAGS) $(SIMBPOBJ) -lm -lz -o $@

matchreadsmain.o: matchreadsmain.cpp matchreadsmain.h UPDATED
	$(CC)  $(INC) $(CFLAGS) -D_BUILD_TIME="\"`cat UPDATED`\"" -c $< -o $@

.cpp.o:  %.o: %.cpp %.h
	$(CC)  $(INC) $(CFLAGS) -c $< -o $@

.c.o:  %.o: %.c %.h
	gcc -O2 -c $< -o $@


clean : 
	rm -fr *.o $(PROGRAMS)
	cd ${SAMTOOLS} && make clean

backup :
	mkdir -p $(BACKUPFOLDER)
	cp *.cpp *.c *.h Makefile *.sh *.pl $(BACKUPFOLDER)
	cp -R ${SAMTOOLS} $(BACKUPFOLDER)
	tar -czf $(BACKUPDIR)/$(STAMP).tgz -C $(BACKUPDIR) $(STAMP)
	rm -fR $(BACKUPFOLDER)
	echo $(BACKUPLOG)
	echo -e "$(STAMP).tgz\t$(comment)" >> $(BACKUPLOG)
