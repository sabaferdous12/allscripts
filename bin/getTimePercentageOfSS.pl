#!/usr/bin/perl -s
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
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If
#   someone else breaks this code, I don't want to be blamed for code
#   that does not work!
#
#   The code may not be sold commercially or included as part of a
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************

use strict;
#use warnings;
use Data::Dumper;
UsageDie() if(defined($::h));

my %resHash;
my $count = 1;

# coil = (~, C, T, S)
# strand = (E, B)
# helix = (H, G, I)

while (my $line=<>) # Reads file from STDIN
{
    chomp($line);
    if ( $line =~ m/^"\S*"/) # Reads line starting with " and non space
    {
        my ($cper, $aper, $sper, $startType);

        my $start = substr ($line, 1, 1); # Takes SS from first frame
        # Checks the type of first frame's SS
        # checks for coil
        if ( ( $start eq "") or ( $start eq "~") or 
                 ( $start eq "C") or ( $start eq "T") or
                     ( $start eq "S") )
           {
               $startType = "coil";
               $cper = getSScount ($startType, $line);
               # A hash storing first frame SS and its
               # percentage over the whole trajectory
               $resHash{$count.$start} = $cper;
           }
        # checks for strand
        elsif ( ( $start eq "E") or ( $start eq "B") )
            {
                $startType = "strand";
                $sper = getSScount ($startType, $line);
                $resHash{$count.$start} = $sper;
                
            }
        # checks for helix
        elsif ( ( $start eq "H") or ( $start eq "G") or
                    ($start eq "I") )
            {
                $startType = "helix";
                $aper = getSScount ($startType, $line);
                $resHash{$count.$start} = $aper;
            }
        $count++;           
    }
}

my ($sum, $countAA) = ( 0, 0);
# printing and average
foreach my $key ( sort {$a <=> $b} keys %resHash)
    {
        $sum = $sum + $resHash{$key};  
        print "$key: $resHash{$key}\n";
        $countAA++;
    }
my $avgTime = $sum/$countAA;
print "Average time for peptide to stay closer to start conformation = $avgTime\n";

# *************************************************************
sub getSScount
{
    my ($start, $simulationSS) = @_;
    my ($coilCount, $strandCount, $alphaCount ) =
        (0, 0, 0);
    my ($cper, $aper, $sper);

    my @ss = split ("", $simulationSS);
    if ( $start eq "coil")
    {
        map { if ( ( $_ eq "~" ) or ( $_ eq "S")or ($_ eq "C")
                       or ($_ eq "T") or ($_ eq "") )
                  {
                      $coilCount++;
                  }
          } @ss;
        # dssp saved frame every 10 ps and its 1000 ns in total
        $cper = (($coilCount*10)/1000000) * 100;
        $cper = sprintf "%.2f", $cper;
        return $cper;
    }
    elsif ( $start eq "strand")
    {
        map{ if ( ( $_ eq "B" ) or ( $_ eq "E") )
                 {
                     $strandCount++;
                 }
         } @ss;
        $sper = (($strandCount*10)/1000000) *100;
        $sper = sprintf "%.2f", $sper;
        return $sper;
    }
    elsif ( $start eq "helix")
    {
        map{ if ( ( $_ eq "H") or ($_ eq "G") or ($_ eq "I") )
                 { 
                     $alphaCount++;
                 }
         } @ss;
        $aper = (($alphaCount*10)/1000000) *100;
        $aper = sprintf "%.2f", $aper;
        return $aper;
    }
    
}


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
Usage:  getTimePercentageOfSS.pl <dsspOutputFile.xpm>

Finds the pecentage of time a residue stays in a parti-
cular SS conformation and returns average over the whole
peptide. It reads an output file from dssp (.xpm) and takes
first frame and then finds its percentage (group of SS) over
the whole trajectory 

__EOF
        exit 0;

    }

    


