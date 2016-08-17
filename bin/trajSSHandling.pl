#!/usr/bin/perl
    #*************************************************************************
    #
    #   Program:    trajSSHandling
    #   File:       trajSSHandling.pl
    #
    #   Version:    V1.0
    #   Date:       01.03.16
    #   Function:   This script extracts frames every 10ps from trajectory and
    #               pass to secondary structure assignment program (pdbsecstr).
    #               It uses SStruc definition for assigning secondry structure.
    #               It then calls bash script trajSS.sh to merge ss column-wise.
    #               And finally it removes pdb files and ss files for each frame.
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
use strict;
#use warnings;
use Getopt::Long qw(GetOptions);

my $trajFile;

    #= $ARGV[0]; # Trajectory file
my $tprFile;

    #= $ARGV[1];  # tpr file
GetOptions
        (
            "t=s" => \$trajFile,
            "s=s" => \$tprFile,
        ) or UsageDie();

my ($name, $ext) = split (/_/, $trajFile);
my $nameLen = length ($name);

# Extract pdb frames from trajectory every 10ps

`gmx trjconv -s $tprFile -f $trajFile -o $name.pdb -dt 10 -sep`;

# Read all PDB files into a directory
use general qw (readDirPDB);
my $dir = ".";

my @PDBfiles = readDirPDB($dir);

# Sort file names numerically - framewise
@PDBfiles = sort{ substr($a, $nameLen) <=> substr($b, $nameLen) } @PDBfiles ;
print join ("\n", @PDBfiles );
my $count=0;

foreach my $pdb ( @PDBfiles)
{
    chomp $pdb;
    print "Calculationg Secondary Structure for frame $count\n";
#    `pdbsecstr -s $pdb | awk '{""; print \$3}' >$name.$count".ss"`;
    `pdbsecstr $pdb | awk '{""; print \$3}' >$name.$count".ss"`;
    
    $count++;
}

# Calling bash script to generate secondary structure of each frame along
# the trajectory

`bash ~/allscript/bin/trajSS.sh $name`;
# Adding " at start and end of each line to make same format as dssp
`perl -p -e 's/^/"/ and s/\$/"/' $name"_temp.dat" >$name"_traj_ss.dat"`;

`rm $name"_temp.dat"`;
`rm *.pdb`;
`rm *.ss`;


#*************************************************************************
# UsageDie()
# ----------
# Prints a usage message and exits
#
# 08.12.15 Original   By: ACRM
sub UsageDie
    {
        print <<__EOF;



trajSSHandling V1.0 (c) 2015, UCL, Saba Ferdous
Usage:  trajSSHandling.pl -t <trajectory file> -s <tpr file>
                
This script extracts frames every 10ps from trajectory and
pass to secondary structure assignment program (pdbsecstr).
It uses SStruc definition for assigning secondry structure.
It then calls bash script trajSS.sh to merge ss column-wise.
And finally it removes pdb files and ss files for each frame
and provides .dat file containing secondary structure of each
residue over whole trajectory.


__EOF
        exit 0;
        
    }
    
