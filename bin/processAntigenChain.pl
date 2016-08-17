#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    processAntigenChain
#   File:       processAntigenChain.pl
#
#   Version:    V1.0
#   Date:       04.01.16
#   Function:   Restoring Antigen chains from PDB and relabeling as 1 if
#               antigen has L or H chain label
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
use general qw (readDirPDB);
use pdb qw (get_pdb_path);
use antigenProcessing qw (getAntigenChains);
use Cwd;
UsageDie() if(defined($::h));

my $dir = getcwd ();

my @dirFiles = readDirPDB($dir);
 
my @antigenPDB;
my $flag = 0;
`mkdir -p processedAntigen`;

foreach my $pdbFile (@dirFiles)
{
    chomp $pdbFile;
    print "Processing $pdbFile ...\n";
    
    my $abFileTemp = "tempAB.pdb";
    my $abFile = "pro^$pdbFile";
    my $antigen = "tempAg.pdb";
    
    # Get 1AFV from 1AFV_1.pdb
    my ($pdbCode, $ext) = split('_',  $pdbFile);

    # Get original PDB path from local PDB
    my $pdbPath = get_pdb_path ($pdbCode);
    
    my @antigenChains = getAntigenChains($pdbFile);
    open (my $IN, '<', $pdbFile) or die "Error - Cannot open input file $!\n ";
    open (my $ABTMP,'>>', $abFileTemp) or die "Error - Cannot open file for writing\n";
    open (my $AB,'>>', $abFile) or die "Error - Cannot open file for writing\n";
    
#    open (my $AG, '>', $antigen) or die "Error....\n";
    
    # Write antibody chains in this file
    while ( my $line = <$IN>) {
        
        if ( ( $line =~ m/^REMARK/) or ( ( substr ($line, 21, 1) eq "L") and ($line =~/^ATOM|^TER/) )
                 or ( ( substr ($line, 21, 1) eq "H") and ($line =~/^ATOM|^TER/) ) ) {
            print {$ABTMP} $line;
        }
    }
    
    #close $ABTMP;
    
    # obtaining orignal antigen chain from local PDB file
    # this is to include any HETATM records in start, middle or end of antigen chain
    my $count = 1;    
    foreach my $ag (@antigenChains)
    {
        open (my $AG, '>', $antigen) or die "Error....\n";
        
        @antigenPDB = `pdbatoms $pdbPath | pdbgetchain $ag`;
       
        
        # If antigen chain has L or H chain label then rename it as A
        if ( ( $ag eq "L" ) or ($ag eq "H") ) {
            map { if ( $_ =~ /\s+$ag\s+/)  {print {$AG} $_;}} @antigenPDB;
            `pdbrenum -c $count -n -d $antigen $count.pdb`;
            # Update header with $count chain label
            $ag = lc($ag);
            
            my $searchStr = "REMARK 950 CHAIN A    $ag";
            
            my $tempStr = "REMARK 950 CHAIN A    $count";
            print $searchStr;
            
            `sed -i -- 's/$searchStr/$tempStr/g' $abFileTemp`;
            `cat $count.pdb >> $abFileTemp`;
            $ag = uc($ag);
            
            @antigenPDB = ();
            `unlink $antigen`;
            `unlink $count.pdb`;
            $count++;
        }
        
        # For normal chain labels
        else {
            # print all lines except master and END record
        #    map { if ($_ !~ /MASTER|END|HOH/)  {print {$ABTMP} $_;}} @antigenPDB;
            open (my $AG2, '>', $ag.".pdb") or die "Error....\n";
            map { if ( $_ =~ /\s+$ag\s*/)  {print {$AG2} $_;}} @antigenPDB;
            `cat $ag.pdb >>$abFileTemp`;
            `unlink $ag.pdb`;
       
        }
                
        @antigenPDB = ();
        `cp $abFileTemp $abFile`;
    }
    
    close $AB;
    close $IN;
    `unlink $abFileTemp`;
    `unlink $antigen`;
#    last;
    
}
#`unlink tempAg.pdb`;

`mv *pro^* processedAntigen`;

# All the new files now have extension pro^:
# To remove that and restore original name, run the following bash command


# `for file in * ; do mv -v "$file" "${file#*^}"; done;` 
# `for file in * ; do sed '/MASTER        0    0    0*/d' "$file" | sed '/"END"/d' >"${file#*^}" ; unlink $file;  done;` 
    
#*************************************************************************
# UsageDie()
# ----------
# Prints a usage message and exits
#
# 08.12.15 Original   By: ACRM
sub UsageDie
{
    print <<__EOF;
processAntigenChain V1.0 (c) 2015, UCL, Saba Ferdous
Usage:  processAntigenChain.pl

processAntigenChain reads all PDB files from current directory
and process antigen chain from antibody-antigen complex by re-
storing it from orignal PDB. It also re-labels chain ID if an
antigen has L or H chain label.
}
__EOF
    exit 0;

}
        
