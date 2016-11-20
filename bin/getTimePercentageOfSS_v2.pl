#!/usr/bin/perl
#*************************************************************************
#
#   Program:    getTimePercentageOfSS
#   File:       getTimePercentageOfSS.pl
#
#   Version:    V1.0
#   Date:       04.01.16
#   Function:   Findings the pecentage of time a residue stays in a parti-
#               cular SS conformation and returns average over the whole
#               peptide. It reads an output file from dssp (.xpm) and takes
#               first frame and then finds its percentage (group of SS) over
#               the whole trajectory
#
#   Copyright:  (c) Saba Ferdous, UCL, 2015
#   Author:     Saba Ferdous
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      saba@bioinf.org.uk
#
#*************************************************************************
#
#   This program is not in the public domain but he code may be modified
#   as required, but any modifications must be
#   documented so that the person responsible can be identified. If
#   someone else breaks this code, I don't want to be blamed for code
#   that does not work!
#
#
#*************************************************************************
#   Usage: ./getTimePercentageOfSS.pl -n 10 -ls 500000 -s 1 -lp 15 -f 2W9EMtPepM4_traj_ss.dat 2 12 13

use strict;
#use warnings;
use Data::Dumper;
use Getopt::Long qw(GetOptions);

my $stepSize;
my $simuLength;
my $inputFile;
my $start;
my $pepLength;

GetOptions
    (
    "n=i" => \$stepSize,
    "ls=i" => \$simuLength,
    "s=i" => \$start,
    "lp=i" => \$pepLength,
    "f=s" => \$inputFile,
) or UsageDie();

open (my $IN, '<', $inputFile);

my %resHash;
my $count = 1;
my @pepSeq;

# coil = (~, -,  C, c,  T, t,  S, s)
# strand = (E, e, B, b)
# helix = (H,h, G, g, I, i)

#################
# Collecting the peptide initial conformation using first frame SS
while (my $line=<$IN>) # Reads file from STDIN
{
    if ( $line =~ m/^"\S*"/ ) # Reads line starting with " and non space
    {
        my $start = substr ($line, 1 , 1); # Takes SS from first frame
        push (@pepSeq, $start);
    }
}

# These are the positions/AA at SS element ends which you want to count
# as stable even if they are loosing the conformation
# e.g -hHHHHHHhTt, I would choose position 2 and 9 to be replaced by X and
# then allow [hHgGtT-~cCsSiI]
print "Peptide Sequence:\n@pepSeq\n";
print "Please Enter positions to be replaced and hit Ctl-d\n";

my @replaceLIST = <STDIN>;

seek $IN, 0, 0;

# Replacing amino acid at chosen positions by letter X
foreach my $replaceList( @replaceLIST)
{
    chomp $replaceList;
    splice @pepSeq, $replaceList-1, 1, "X";
}
print "The new peptide sequence:\n@pepSeq\n";
################


my $count = 0;
#my $count = $start-1; 
my $ssCount=0;
my $sum = 0;
my $ssper;
 
while (my $line2=<$IN>)
{
    chomp($line2);
    my $first = $pepSeq[$count];
    $count++;

    my @ss = split ("", $line2);
    if ( $line2 =~ m/^"\S*"/ ) # Reads line starting with " and non space
    {
        $ssCount=0;
        if ( ( $first eq "") or ( $first eq "~") or ( $first eq "-") or
                 ( $first eq "C") or ( $first eq "c") or
                     ( $first eq "S") or ( $first eq "s" ) or
                         ( $first eq "T") or ( $first eq "t" ) )
            {
            map { if ( ( $_ eq "~" )  or ( $_ eq "-")or ($_ eq "")
                         or ($_ eq "S") or ($_ eq "s") or
                             ( $_ eq "C") or ( $_ eq "c") or
                             ($_ eq "T") or ($_ eq "t") )
                      {
                          $ssCount++;
                      }
              } @ss;
                
        }#########


        if ( $first eq "X")
        {
            map { if ( ( $_ eq "~" )  or ( $_ eq "-")or ($_ eq "") or 
                           ($_ eq "S") or ($_ eq "s") or
                               ( $_ eq "C") or ( $_ eq "c") or
                                   ($_ eq "T") or ($_ eq "t") or
                                       ($_ eq "h") or ($_ eq "H") or
                                           ($_ eq "h") or ($_ eq "G") or
                                               ($_ eq "i") or ($_ eq "I") )
                      {
                          $ssCount++;
                      }
              } @ss;
        }#####

        if ( ($first eq "h") or ($first eq "H") or
                 ($first eq "g") or ($first eq "G") or
                     ($first eq "i") or ($first eq "I") ) {
            
            map { if ( ($_ eq "h") or ($_ eq "H")or
                           ($_ eq "g") or ($_ eq "G") or
                               ($_ eq "i") or ($_ eq "I") )
                      {
                          $ssCount++;
                      }
              } @ss;
        }######


        if ( ($first eq "e") or ($first eq "E") or
                     ($first eq "b") or ($first eq "B") ) {
            
            map { if ( ($_ eq "e") or ($_ eq "E") or
                           ($_ eq "b") or ($_ eq "B") )
                      {
                          $ssCount++;
                      }
              } @ss;
        }###### 
        
        
    }
    
    
### Percentage of each residue
    # dssp saved frame every 10 ps and its 1000 ns in total
    # $stepSize = 10ps, $simuLength = 1000000ps 
    $ssper = (($ssCount*$stepSize)/$simuLength) * 100;
    $ssper = sprintf "%.2f", $ssper;
    $resHash{$count.$first} = $ssper;
    #    $sum = $sum+$ssper;
    #    print "$first: $ssper\n";
}


my @perVal;
my @perKey;
my $avgTime; 
foreach my $key ( sort {$a <=> $b} keys %resHash)
{
    push (@perVal, $resHash{$key} );
    push (@perKey, $key );
  #  print "$key: $resHash{$key}\n";              
}

my $startLoop = $start-1;
my $endLoop = $pepLength + $startLoop;

my $countAA=0;
for ( my $i = $startLoop; $i < $endLoop; $i++) {
    $sum = $sum + $perVal[$i];
    $countAA++;
    print "$perKey[$i]: $perVal[$i]\n";
}

#print "SUM=$sum\n";

$avgTime= sprintf "%.2f", $sum/$countAA,"\n";
print "Average time for peptide to stay closer to start conformation = $avgTime\n";



#*************************************************************************
# UsageDie()
# ----------
# Prints a usage message and exits
#
# 08.12.15 Original   By: ACRM
sub UsageDie
    {
        print <<__EOF;


getTimePercentageOfSS V1.0 (c) 2015, UCL, Saba Ferdous
Usage:  getTimePercentageOfSS.pl -n <step size in ps> -ls <simulation length in ps>
        -s <start position> -lp <length of peptide>
        -f <dsspOutputFile.xpm>

Finds the pecentage of time a residue stays in a parti-
cular SS conformation and returns average over the whole
peptide. It reads an output file from dssp (.xpm) and takes
first frame and then finds its percentage (group of SS) over
the whole trajectory 

__EOF
        exit 0;

    }

    


