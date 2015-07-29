#!/usr/bin/perl
#*************************************************************************
#
#   Program:    constructNonRedundantData
#   File:       constructNonRedundantData.pl
#   
#   Version:    V1.1
#   Date:       22.04.14
#   Function:   Take first element from file containing redundanat clusters  
#               and construct non-redundant dataset 
#   Usage:      ./constructNonRedundantData.pl <Redundant cluster file> 
#                  <output directory with full path>
#
#   Copyright:  (c) UCL, Saba Ferdous, 2014
#   Author:     Miss Saba Ferdous
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      saba@bioinf.org.uk
#               
#*************************************************************************
use strict;
use warnings;
use File::Copy;

my $infile = $ARGV[0]; # List of redundant clusters
my $outputDir = $ARGV[1]; # Directory to store non-redundant data

open(my $IN, $infile) or die "Can not open $!\n";
my @nrFiles;
while (<$IN>)
{
    my ($firstElemInCluster) = split(',', $_);
    chomp $firstElemInCluster; 
 
    my $firstElemPDBfile =  $firstElemInCluster.".pdb";
    push (@nrFiles, $firstElemPDBfile);
}

# Copying non-redundant data into given directory
foreach my $nrFile (@nrFiles)
{
    chomp $nrFile;
    my $dest = "../$outputDir";
    if (-e $nrFile)
    {
        copy ($nrFile, $dest);
    }
    
    else
    {
        print "$nrFile does not exists\n"; 
    }
}
