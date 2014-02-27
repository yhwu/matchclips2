## Update

This is an updated version of [matchclips](https://github.com/yhwu/matchclips) to precisely locate the break points of structure variations. This update mainly focuses on speed and reproducibility, and is intended for real data large files. A 70x coverage whole genome sequencing BAM file should take about 5-8 hours on a single cpu and use about 2GB memory. It is not really helpful to use more than 4 cores. Performance suffers if the computing node is disk IO busy. 

This program works as follows:

1. use read depths to detect longer CNVs and test them ( to be added later, [rsicnv](https://github.com/yhwu/rsicnv)! );
2. use paired end distances to detect longer structure variations and test them;
3. use split reads to detect shorter structure variations and test them.

The tests include:

1. read depth test --- a CNV should have a read depth either lower than 2/3 or higher than 5/4 of its neighbor's. These two ratios are calculated assuming two chromosomes are independent and their fragmentation follows the same normal distribution.
2. read pair test --- there must be enough abnormal pairs properly envelope the break points of a structure variation. Minimum required is 1/3 of expected number of pairs(templates). No reasoning here, purely empirical, just to take mapping difficulty into consideration.
3. matching reads test --- reads on left side break point must overlap with reads on right side break point. Minimum number required is 1/4 of expected. No reasoning here, purely empirical, with mapping difficulty and length of overlap taken into consideration.

Please cite MATCHCLIPS by ```doi: 10.3389/fgene.2013.00157```.

## Output fields
```
  1. CHROM
  2. BEGIN, lower position on reference, 1 based
  3. END,  higher position on reference, 1 based
  4. DEL or DUP
     The break points are either DEL or DUP, which does not necessarily
     mean deletion or duplication. For example, a translocation will give
     both types of break points, but no copy number change. 
     For DEL, the 5' break point is BEGIN and 3' is END,
     For DUP, the 5' break point is END and 3' is BEGIN.
     note: matchclips does not make calls on what type of structure 
     variation those break points or combinations of break points represent.
  5. LENGTH = 3' break point - 5' break point
  6. UN:n
     Sometimes, moving both break points simultaneously will not change
     the variation sequence due to repeats in the reference. This number
     tells how many bases are repeated at the two break points. please
     refer to the paper or the old matchclips for clarification.
  7. RD:n1;n2;n3:s read depth information
     n1: average read depth within (BEGIN-dx, BEGIN)
     n2: average read depth within (END, END+dx)
     n3: average read depth within (BEGIN, END)
      s: score of read depth
          0, not likely,
          1, kind of, 
          2, likely,
          3, very likely,
          4, maybe double deletion.
     note: dx=2*(END-NEGIN), l_qseq < dx < 4*l_qseq
  8. RP:n1;n2;n3:s read pair information
     n1: normal pairs across BEGIN
     n2: normal pairs across END
     n3: pairs that properly envelope 5' and 3' break points
      s: strength of signal
          0, not likely,
          1, kind of,
          2, likely
          3, very likely,
     note: template of length within mean +- 5 s.d. is considered normal, 
     please read on for clarification.
  9. MR:n1;n2;n3:s matching reads information
     n1: reads on BEGIN side that match reads on END side
     n2: reads on END side that match reads on BEGIN side
     n3: average edit distance between the merged reads and the corresponding
         reference,
      s: strength of signal
          0, not likely,
          1, kind of
          2, likely,
          3, very likely,
     note: reads matching is done in the direction from 5' to 3'.
 10. SR:n1;n2 matching softclipped reads information
     n1: averaged edit distance
     n2: number of read pairs that match, e.g., if 3 reads on 5' side match 4
         reads on 3' side, this number will be 12.
     note: the reads have long softclipped part(S>=10) 
 11. Q0:n1;n2 percentage of low mapping quality within BEGIN and END
     n1: percentage of mapq==0
     n2: percentage of mapq<=10
```
Note: 
 Any negative number means the field is not calculated, either because it is not necessary, or not useful.
 Along the main ```$output``` file, ```$output.weak``` contains CNVs of weaker signal.
 Short indels(<10 bases) are detected too. But, seriously, short indels ought to be detected directly by mapping software.     

## Usage
```
[yhwu@penn matchclip2]$ ./matchclips 
Usage:
  matchclips <options> -f REFFILE -b BAMFILE [REGIONS]

Options:
  -t  INT  number of threads, INT=1 
  -e  INT  max allowed mismatches when matching strings, INT=2 
  -l  INT  minimum length of overlap, INT=25 
  -s  INT  minimum number of soft clipped bases, INT=10 
  -q  INT  minimum mapping score, INT=10 
  -Q  INT  minimum base read quality, INT=0 
  -L  INT  check reads before and after INT bases, INT=auto 
  -L  0    pair end mode only 
  -se      single end mode, do not use pair end distances
  -pe INT INT provide insert and s.d. of insert, otherwise calculate them
  -o  STR  outputfile, STR=STDOUT 
   REGION  if given should be in samtools's region format
```

## Code example
```
./matchclips                                                 #help
./matchclips -f hg19.fasta -b A.bam -o A.txt 
./matchclips -t 4 -f hg19.fasta -b A.bam chr1 chr2 -o A.txt  #use 4 cores and process chr1 and chr2
./matchclips -se -L 1000000 -f hg19.fasta -b A.bam -o A.txt  #single end mode and check reads matching within 1000000
                                                             #this is equivalent to the original matchclips
./matchclips -L 0 -f hg19.fasta -b A.bam -o A.txt            #paired end mode only
```

## Installation
-. github way:
```
git clone https://github.com/yhwu/matchclips2.git
cd matchclips2 
make
```
-. manual download:
```
wget https://github.com/yhwu/matchclips2/archive/master.zip -O matchclips2.zip
unzip matchclips2.zip
cd matchclips2-master
make
```

## On target sequencing, tumor, exom

Sometimes we like to do target sequencing on tumor samples because there must be something wrong, and we end up with millions of reads within a small region but there is no signal. There might be several reasons. 
 
 First, in the tumor sample, there might be only a small fraction of cells that are really bad, so signals are rare to begin with. 

 Second, the enrichment process is selective and it favors good DNAs because that is what we know as the reference, hg19 or whatever. The baits bond to matched DNA sequences but not the altered ones. If baiting efficiency is 1 for good DNA, the baiting efficiency for varied sequence templates is roughly ```1 - 2*baitLength/templateLength```. So, enriching the target many times effectively subdues the bad DNA. 

 Thirdly, the bad cells may have DNAs damaged to different degrees and the results are smeared and clean signals may not be possible. Yet, all methods are based on some assumptions and usually require clean signals. The varied pieces of templates are also difficult to map, making it difficult for paired end and split read methods. If the variation is long, e.g., a long deletion or duplication, read depth methods can still be used. This is why we must add read depth method into the package. However, it is difficult because enrichment messes up every thing.

For the reasons above, we do not throw away weak signal. They are saved in ```$output.weak``` along with the main output. Your enriched target sequencing may only contain trace evidence.

## On read depth 

It is a misconception to think read depth information is only useful for longer CNVs. For example, if we have a deletion somewhere, we ought to see a sudden drop of depth at the 5' break point and a sudden jump of depth back to normal at the 3' break point. Within a shorter distance, the read depth is more likely to stay roughly the same than over a longer distance. So, maybe, read depth information is more reliable for shorter CNVs. For the same reason, I do not think adjusting for GC content is helpful.

## On paired end sequencing

-. paired end method can be precise too

The commonly used Illumina paired end sequencing sequences both segments of a template. A template has two strands, and they are reverse complement to each other. Since the sequencers only read from 5' to 3', so one read(segment) must be on the FORWARD strand and the other must be on the REVERSE COMPLEMENT strand. It is the aligners(bwa, botie, novoalign, etc.)'s job to figure out which is which, and the orientation of a read is coded in the FLAG field in the SAM format. There are four coordinates that are import to describe a pair(template),<br> 
  F1, start of the forward read,<br>
  F2, end of the forward read,<br>
  R1, the start of the reverse complement read,<br>
  R2, the end of the reverse complement read,<br>
and all the four(three) numbers can be extracted from a FORWARD(REVERSE) read, as shown in the figure below.

![The pair end algorithm](https://dl.dropboxusercontent.com/u/39236968/pairendtemplate.svg)

A structure variation concatenates two bases(break points) that are not next to each other on the reference. Templates sampled around the variation still have actual lengths that follow normal distribution just like other templates on normal regions do. But the template lengths calculated from the coordinates will have a different mean, as shown in the figure below. The template lengths given in SAM format are R2-F1 for the FORWARD reads and -(R2-F1) for the reverse reads, and abnormal template lengths can be used to infer structure variations. It is also important to note that, if mapping is done correctly, the F2's will be always lower than the 5' side of the break point, and the R1's will always be higher than the 3' side of the break point. The reads may actually cover the break points, but the parts that pass the break points are usually and should be marked as soft clipped, and those parts don't count in the calculation of the start and the end of a read. See the figure below. So, if there are many templates sampled across the variation, we should be able to obtain the exact break points. This is indeed the case from our experience in developing the program.

![The pair end algorithm](https://dl.dropboxusercontent.com/u/39236968/pairendbreakpoint.svg)
	
-. our paired end test is to "properly" "envelope"
	
Based on the above analysis, our pair end test looks for pairs whose F1's are lower than(to the left on coordinate of) the 5' side break point, and R2's  higher than the 3' side break point, and whose actual template lengths as calculated in the figure fall within 5 s.d. of the mean of normal pairs. So, "properly" means the pairs are rightly sized, and "envelope" means the start and end of a pair are outside of the break points. For example, for a deletion, the start(F1) of a pair is on the lower than the lower break point, and the end(R2) of the pair is higher than the higher break point; for a duplication, the start(F1) of a pair is on the lower than the higher break point, and the end(R2) of the pair is higher than the lower break point. 

## Matching reads test

This is similar to the [matchclips](https://github.com/yhwu/matchclips) algorithm except that we do not require the reads are soft clipped. Given a pair of break points, we collect a set of reads that cover the 5' break point, and a set for the 3' break point. We then, test how many of reads in the left(5') set overlap with how many of the reads in the right set. Since reads that span the break points are hard to map, we require 1/4 of the reads overlap. Two reads are considered overlaped if they overlap at least 35 bases to the ends.


## AUTHOR
Yinghua Wu, Hongzhe Li
Department of Biostistics and Epidemiology, University of Pennsylvania, Philadelphia, PA 19104 

## License
GNU General Public License, version 3.0 (GPLv3)
