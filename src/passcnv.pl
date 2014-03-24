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

my $cnvfile="";
$cnvfile=shift @ARGV if ( $numArgs>=1 ) || exit ;

my ( $chr, $beg, $end, $type, $length, $un,
     $rdL, $rdR, $rd, $rdscore,
     $rpL, $rpR, $rp, $rpscore,
     $mrL, $mrR, $mrE, $mrscore,
     $srE, $sr,
     $q0, $q10 );

open (FIN, $cnvfile ) or die "$cnvfile not found $!\n";
while (<FIN>) {
    s/[\n\r+]$//;
    if (/^#/) {print "$_\n"; next;}
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
    
##### filtering goes here #########
    my $is_false=0;
    
    if ( $rdscore>=0  && $rdscore<=1 ) { $is_false=1; }
    if ( $rpscore>=0  && $rpscore<=1 ) { $is_false=1; }
    if ( $mrscore>=0  && $mrscore<=1 ) { $is_false=1; }
    
    print "$_\n" if ( ! $is_false ) ;
    
##### filtering ends here #########
}

close(FIN);
exit;

sub usage {
    die(qq/
This perl script is used to pass the output result from matchclips so that
users can select variations according to their own creterions. User's code
should be added between the two comment lines. 

##### filtering goes here #########

##### filtering ends here #########

Available variables include

\$chr, \$beg, \$end, \$type, \$length, \$un,
\$rdL, \$rdR, \$rd, \$rdscore,
\$rpL, \$rpR, \$rp, \$rpscore,
\$mrL, \$mrR, \$mrE, \$mrscore,
\$srE, \$sr,
\$q0, \$q10 

Examples: 
          passcnv.pl cnvlist1.txt
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

