#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    findPeptideAntigen
#   File:       findPeptideAntigen.pl
#
#   Version:    V1.0
#   Date:       04.01.16
#   Function:   Identification of PDB files with peptide antigens
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
use warnings;
use List::MoreUtils qw(uniq);
use general qw (readDirPDB);
use antigenProcessing qw (getAntigenChains);

use SFPerlVars;
use Cwd;

UsageDie() if(defined($::h));
              
my $dir = getcwd ();

my @dirFiles = readDirPDB($dir);
my @PDBFiles;

foreach my $pdbFile (@dirFiles)
{
    chomp $pdbFile;
    # Obtain antigen chain labels 
    my @antigenChains = getAntigenChains($pdbFile);

    # Check for peptide chain type
    foreach my $ag ( @antigenChains ) {
        my $chaintype = `chaintype -c $ag -p $pdbFile\n`;
        if ( $chaintype =~ /PEPTIDE/) {
            push (@PDBFiles, $pdbFile);
                       
        }
    }
    @PDBFiles = uniq @PDBFiles;
      
}
print join ("\n", @PDBFiles );

    
#*************************************************************************
# UsageDie()
# ----------
# Prints a usage message and exits
#
# 08.12.15 Original   By: ACRM
sub UsageDie
{
    print <<__EOF;

findPeptideAntigen V1.0 (c) 2015, UCL, Saba Ferdous
Usage:  findPeptideAntigen.pl 
        
findPeptideAntigen reads all PDB files from current directory
and provides a list of PDB files with peptide antigen.

__EOF
    exit 0;
}
