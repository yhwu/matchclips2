PROGRAMS = matchclips tagcnv
SAMTOOLS = samtools-0.1.19
all: $(PROGRAMS)
TITLE = matchclip2
CC =g++
CFLAGS =-g -Wall -O2 
#CFLAGS = -g -Wall -L/usr/lib64
INC=-I${SAMTOOLS}
LIBS =   -L/usr/lib64 -L${SAMTOOLS} -lbam -pthread -lm -lz
BACKUPDIR = ~/projects/$(strip $(TITLE))
BACKUPLOG = $(strip $(BACKUPDIR)/$(TITLE)).log
STAMP = $(strip $(shell  date +'%Y.%m.%d-%H.%M.%S'))
BACKUPFOLDER = $(BACKUPDIR)/$(STAMP)	

TAGCXX =  tagcnvmain.cpp tagcnv.cpp samfunctions.cpp readref.cpp functions.cpp
TAGHDR = $(TAGCXX:.cpp=.h)	
TAGOBJ = $(TAGCXX:.cpp=.o)	
tagcnv : $(TAGOBJ) $(TAGCXX) $(TAGHDR) Makefile
	$(CC) $(CFLAGS) $(TAGOBJ) $(INC) $(LIBS)  -o $@

MATCHCXX =  matchreadsmain.cpp matchreads.cpp preprocess.cpp exhaustive.cpp pairguide.cpp  samfunctions.cpp readref.cpp functions.cpp
MATCHHDR = $(MATCHCXX:.cpp=.h)	
MATCHOBJ = $(MATCHCXX:.cpp=.o)	
matchclips : $(MATCHOBJ) $(MATCHCXX) $(MATCHHDR) Makefile ./${SAMTOOLS}/libbam.a
	$(CC) $(CFLAGS) $(MATCHOBJ) $(INC) $(LIBS) -o $@

time:
	date -u "+%a %b %d %H:%M:%S %Y" > UPDATED

./${SAMTOOLS}/libbam.a:
	cd ${SAMTOOLS} && make libbam.a

matchreadsmain.o: matchreadsmain.cpp matchreadsmain.h UPDATED
	$(CC) -c $(CFLAGS) $(INC) -D_BUILD_TIME="\"`cat UPDATED`\"" $< -o $@

.cpp.o:  %.o: %.cpp %.h
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

.c.o:  %.o: %.c %.h
	gcc -c -wall -O2 $< -o $@


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
