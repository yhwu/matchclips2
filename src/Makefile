PROGRAMS = matchclips
SAMTOOLS = samtools-0.1.19
all: $(PROGRAMS)
TITLE = matchclip2
CC = g++
CFLAGS = -g -O2 -Wall -L/usr/lib64
#CFLAGS = -g -Wall -L/usr/lib64
INC=-I${SAMTOOLS}
LIBS =  -L${SAMTOOLS} -lbam -pthread -lm -lz
BACKUPDIR = ~/projects/$(strip $(TITLE))
BACKUPLOG = $(strip $(BACKUPDIR)/$(TITLE)).log
STAMP = $(strip $(shell  date +'%Y.%m.%d-%H.%M.%S'))
BACKUPFOLDER = $(BACKUPDIR)/$(STAMP)	

MATCHCXX =  matchreadsmain.cpp matchreads.cpp preprocess.cpp exhaustive.cpp pairguide.cpp  samfunctions.cpp readref.cpp functions.cpp
MATCHHDR = $(MATCHCXX:.cpp=.h)	
MATCHOBJ = $(MATCHCXX:.cpp=.o)	
matchclips : $(MATCHOBJ) $(MATCHCXX) $(MATCHHDR) Makefile ./${SAMTOOLS}/libbam.a
	$(CC) $(CFLAGS) $(MATCHOBJ) $(INC) $(LIBS) -o $@

./${SAMTOOLS}/libbam.a :
	cd ${SAMTOOLS} && make libbam.a

.cpp.o:  %.o: %.cpp %.h
	$(CC)  $(INC) $(CFLAGS) -c $< -o $@

.c.o:  %.o: %.c %.h
	gcc -O2 -c $< -o $@


clean : 
	rm -fr *.o $(PROGRAMS)
	cd ${SAMTOOLS} && make clean

backup :
	mkdir -p $(BACKUPFOLDER)
	cp *.cpp *.h Makefile $(BACKUPFOLDER)
	cp -R ${SAMTOOLS} $(BACKUPFOLDER)
	tar -czf $(BACKUPDIR)/$(STAMP).tgz -C $(BACKUPDIR) $(STAMP)
	rm -fR $(BACKUPFOLDER)
	echo $(BACKUPLOG)
	echo -e "$(STAMP).tgz\t$(comment)" >> $(BACKUPLOG)
