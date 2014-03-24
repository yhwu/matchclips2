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

usage() if ( $numArgs<1 ) ;

my $precise=0;
my $minLen=0;
my $precision=22;
for( my $i=0; $i<@ARGV; $i++) {
    if ( $ARGV[$i] eq "-p" ) {
	$precise=1;
	shift @ARGV;
	$precision=shift @ARGV;
    }
    if ( $ARGV[$i] eq "-l" ) {
	shift @ARGV;
	$minLen=shift @ARGV;
    }
}
#print $precision,"\t",$minLen,"\n"; exit(0);

my $fastq=$ARGV[0];
my $cnvfile="";
$cnvfile=$ARGV[1] if ( $#ARGV > 0 );

print STDERR "#reads file : ",$fastq,"\n";
print STDERR "#cnv   file : ",$cnvfile,"\n";
#exit;

my ($prechr, $prep1, $prep2)=(0,1,2);

my (@chr, @p1, @p2, @pos5, @pos3, %idx, @is_detected);

open (FIN, $fastq ) or die "$fastq not found $!\n";
my $linein=<FIN>;
my $is_my_sim_fq = $linein =~ m/^@/ ? 1:0;
seek(FIN, 0, 0) or die "Can't seek to beginning of file: $!";
$.--;
while ( $linein=<FIN>) {
    my @cell;
    if ( $is_my_sim_fq ) {
	next if ( $.%4 != 1 );
	@cell = split (/_/,$linein);
    }
    else { @cell = split (/\s/,$linein); }
    
    next if $cell[2]==0; 
    
    $cell[0] =~ s/@//;
    $cell[0] =~ s/chr//;
    $cell[0] = 23 if $cell[0] eq "X" ;
    $cell[0] = 24 if $cell[0] eq "Y" ;
    ($cell[1], $cell[2])=($cell[2], $cell[1]) if ( $cell[1] > $cell[2] );
    next if ( $cell[2]-$cell[1] < $minLen );
    next if ( $cell[0] eq $prechr && $cell[1] eq $prep1 && $cell[2] eq $prep2 ) ;
    $prechr=$cell[0];
    $prep1=$cell[1];
    $prep2=$cell[2];
    
    push @chr,  $cell[0];
    push @p1,  $cell[1];
    push @p2,  $cell[2];
    push @pos5,  $cell[0]*1E10+$cell[1];
    push @pos3,  $cell[0]*1E10+$cell[2];
    push @is_detected, 0;
    $idx{$pos5[-1]} = $#pos5;
#    print $cell[0],"\t",$cell[1],"\t",$cell[2],"\t", $pos5[-1],"\n";
}
close(FIN);
print STDERR "-------------------------\t",$#pos5+1,"\n";

exit if ( $cnvfile eq "" );

open (FIN, $cnvfile ) or die "$cnvfile not found $!\n";
my $falsepossitive=0;
my $myfound=0;
my $ismatchclips=0;
while ( $linein=<FIN>) {
    $linein =~ s/[\n\r+]$//;
    $ismatchclips=1 if ( $linein =~ m/UNCERTAINTY/ );
    if ( $linein =~ /^\#/ ) {
	print $linein,"\tOVERLAP","\n";
	next;
    };
##    next if ( $linein !~ /BP_/ );

    my @cell = split (/\s/,$linein);
    
    $cell[0] =~ s/chr//;
    $cell[0] = 23 if $cell[0] eq "X" ;
    $cell[0] = 24 if $cell[0] eq "Y" ;
    next if ( $cell[0] !~ (/^\d+$/) );

    ($cell[1], $cell[2])=($cell[2], $cell[1]) if ( $cell[1] > $cell[2] );

    next if ( $cell[2]-$cell[1] < $minLen );
    my $cnvlen=$cell[2]-$cell[1]+1;

#    next if ( $cell[2]-$cell[1] < 5 );

    my $pos=$cell[0]*1E10+$cell[1];
    my $pos3end=$cell[0]*1E10+$cell[2];
    my $dx=0;
    $dx=$cell[4] if $ismatchclips ;
    $dx=0;
#    if ( $#cell>11){
#	$dx=$cell[7]-$cell[8]-$cell[11];
#	$dx=$cell[4] if $cell[4]>$dx;
#	$dx=$cell[4];
#    } 
    
    my $is_there=0;
    for( my $i=0; $i<@pos5; $i++) {
	if ( $precise ) {
	    if ( abs($pos-$pos5[$i])<=$precision+$dx && 
		 abs($pos3end-$pos3[$i])<=$precision+$dx ) {
		$is_there=1;
		$is_detected[$i]=1;
	    }
	}
	else {
	    my $max5end= $pos>$pos5[$i] ? $pos:$pos5[$i];
	    my $min3end= $pos3end<$pos3[$i] ? $pos3end:$pos3[$i];
	    if ( $max5end-$dx<$min3end+$dx ) {
		$is_there=1;
		$is_detected[$i]=1;
	    }
	}
    }
    
    $falsepossitive+=1-$is_there if $cnvlen>5; #simulation introduces short indels
#    $falsepossitive+=1-$is_there;
    $myfound+=$is_there;
    print $linein,"\t",$is_there,"\n";
}
my $num_found=0;
$num_found+=$_ for( @is_detected );

#print "Found $num_found out of ",$#pos5+1,"\t";
#print "False possitive ",$falsepossitive,"\n";
my $filename=$cnvfile;
print STDERR "\n#filename\tprecision\tnum_found\tfalsepositive\trefset\tFound\n";
print STDERR $filename,"\t",$precise*$precision,"\t",$num_found,"\t",$falsepossitive,"\t",$#pos5+1,"\tFound\n";
close(FIN);

exit;

sub usage {
    die(qq/
Examples: 
          checkbp.pl cnvlist1.txt cnvlist2.txt       #overlap
          checkbp.pl -p 3 cnvlist1.txt cnvlist2.txt  #precise to 3 bases
\n/);
}
