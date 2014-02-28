#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

my @bases=split(/\//, $0);
my $me=$bases[-1];
my $numArgs = $#ARGV + 1;
my @ARGV2=@ARGV;

usage() if ( $numArgs<2 ) ;

my($cnvfile, $annofile)= ($ARGV[0], $ARGV[1]); 
print STDERR "#cnv   file : ",$cnvfile,"\n";
print STDERR "#anno  file : ",$annofile,"\n";

my (@chr, @p1, @p2, @anno);

open (FIN, $annofile ) || die "$annofile not found $!\n";
my $anno_comment="";
my $linein=<FIN>;
if ( $linein =~ /Chr,Start,End/ ) {
    $linein =~ s/[\n\r+]$//;
    $anno_comment=( split /,/, $linein, 6 )[5] ;
}
print STDERR "#anno  field : ",$anno_comment,"\n";

while ( $linein=<FIN>) {
    $linein =~ s/[\n\r+]$//;
    next if $linein =~ /^#/; 
    my( $rname, $start, $end, $ref, $alt, $txt) = ( split /,/, $linein, 6 ) ;
    push @chr,  $rname;
    push @p1,  $start;
    push @p2,  $end;
    push @anno,  $txt;
}
close(FIN);
print STDERR "#anno  field : ",$#anno+1,"\n";

open (FIN, $cnvfile ) || die "$cnvfile not found $!\n";
my $idx=0;
while ( $linein=<FIN>) {
    $linein =~ s/[\n\r+]$//;
    if ( $linein =~ /^#/ ) {
	print $linein,"\t",$anno_comment,"\n";
	next;
    }
    my( $c, $beg, $end ) = ( split /\s/, $linein, 4 ) ;
    my $found=-1;
    for( my $i=$idx; $i<@anno; ++$i) {
	$found=$i if ( "$c|$beg|$end" eq "$chr[$i]|$p1[$i]|$p2[$i]" ); 
	last if ( $c eq $chr[$i] && $beg < $p1[$i] );
    }
    if ( $found>0 ) {
	$idx=$found;
	print $linein,"\t", $anno[$found],"\n";
    }
    else { print $linein,"\t#\n"; }
}
close(FIN);


exit;

sub usage {
    die(qq/
Examples: 
          append_anno.pl cnvfile.txt annofile.txt
\n/);
}
