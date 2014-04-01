#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

my @bases=split(/\//, $0);
my $me=$bases[-1];
my $numArgs = $#ARGV + 1;

usage() if ( $numArgs<1 ) ;

my ($bam, $cnvfile)=("", "");

for( my $i=0; $i<@ARGV; $i++) {
    if ( $ARGV[$i] eq "-b" ) { shift @ARGV; $bam=shift @ARGV; }
    if ( $ARGV[$i] eq "-cnv" ) { shift @ARGV; $cnvfile=shift @ARGV; }
}

usage() if ( $bam eq "" || $cnvfile eq "" ) ;

my ( $chr, $beg, $end, $type, $length, $un,
     $rdL, $rdR, $rd, $rdscore,
     $rpL, $rpR, $rp, $rpscore,
     $mrL, $mrR, $mrE, $mrscore,
     $srE, $sr,
     $q0, $q10 );

open (FIN, $cnvfile ) or die "$cnvfile not found $!\n";
print "new\ngenome hg19\n";
print "load $bam\n";
print "snapshotDirectory ./snapshots\n\n";

while (<FIN>) {
    s/[\n\r+]$//;
    next if (/^#/);
    next if ( length($_)<3 );
    
    ( $chr, $beg, $end, $type, $length, $un,
      $rdL, $rdR, $rd, $rdscore,
      $rpL, $rpR, $rp, $rpscore,
      $mrL, $mrR, $mrE, $mrscore,
      $srE, $sr,
      $q0, $q10 ) = 
	  $_ =~ /^(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\tUN:(.*?)\tRD:(.*?);(.*?);(.*?):(.*?)\tRP:(.*?);(.*?);(.*?):(.*?)\tMR:(.*?);(.*?);(.*?):(.*?)\tSR:(.*?);(.*?)\tQ0:(.*?);(.*?)$/;
    if ( cnv_format() ne $_ ) {  ## double check
	print STDERR "input not passed correctly\n";
	print STDERR $_,"\n";
	print STDERR cnv_format(),"\n";
	exit;
    }
    
    printigv() ;
    next;

##### filtering goes here #########
    my $is_false=0;
    
    if ( $rdscore>=0  && $rdscore<=1 ) { $is_false=1; }
    if ( $rpscore>=0  && $rpscore<=1 ) { $is_false=1; }
    if ( $mrscore>=0  && $mrscore<=1 ) { $is_false=1; }
    
    printigv() if  ( ! $is_false ) ;
#    print "$_\n" if ( ! $is_false ) ;
    
##### filtering ends here #########
}

close(FIN);

print "exit\n";

exit;

sub usage {
    die(qq/
This perl script generates igv batch script to make igv snapshots around
CNVs. By default, regions twice as long as the CNV is shown before and after 
the CNV. 

Examples: 
          igvbatch.pl -b input.bam -cnv cnvlist1.txt > batch.txt
          igv.sh -g hg19 -b batch.txt
\n/);
}

sub cnv_format {
    return 
	"$chr\t$beg\t$end\t$type\t$length\tUN:$un\t".
	"RD:$rdL;$rdR;$rd:$rdscore\t".
	"RP:$rpL;$rpR;$rp:$rpscore\t".
	"MR:$mrL;$mrR;$mrE:$mrscore\t".
	"SR:$srE;$sr\t".
	"Q0:$q0;$q10";
}

sub printcnv {
    print "$chr\t$beg\t$end\t$type\t$length\t$un\t";
    print "RD:$rdL;$rdR;$rd:$rdscore\t";
    print "RP:$rpL;$rpR;$rp:$rpscore\t";
    print "MR:$mrL;$mrR;$mrR:$mrscore\t";
    print "SR:$srE;$sr\t";
    print "Q0:$q0;$q10\n";
}

sub printigv {
    
    my $ichr=$chr;
    my $ibeg=$beg;
    my $iend=$end;
    swap($ibeg, $iend) if $ibeg > $iend;
    $ichr="chr$ichr" if ($ichr !~ /\D/); 
    my $dx = $iend - $ibeg;
    $dx=100 if ( $dx<100 );
    $dx=5000 if ( $dx>5000 );
    my $left=$ibeg-$dx;
    my $right=$iend+$dx;
    $dx = $iend - $ibeg;
    my $snapshot="$chr\_$ibeg-$iend-$type-$dx.png";
    print "region $ichr $ibeg $iend\n";
    print "goto $ichr:$left-$right\n";
    # print "sort base\n";
    print "collapse\n";
    print "snapshot $snapshot\n";
    print "\n";
}
