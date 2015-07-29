#!/usr/bin/perl
#*************************************************************************
#
#   Program:    getRedundantAntibodyClusters
#   File:       getRedundantAntibodyClusters.pl
#   
#   Version:    V1.1
#   Date:       22.04.15
#   Function:   For the given directory for antibody PDB structures/files, 
#               computes the clusters of redundant antibodoies
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
#
#   Description:
#   ============
#   The script finds the redudant clusters of the antibodies with identical 
#   sequences. It works by taking substring of protein (Middle part) being 
#   compared and then compares ends for 100% identity
#
#*************************************************************************
#
#   Usage:
#   ======
#  ./getRedundantAntibodyClusters.pl 
#   (Script should be run from the desired directory)
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.1   22.04.15 Original
#
#*************************************************************************

use strict;
#use warnings;
use Data::Dumper;
use Cwd;
use File::Basename;

use redund qw(getAntibodySeq aa3to1 readDir splitString
getEndsSeq compareStarts compareEnds checkSubstring printClusters 
isIdenticalEnds);

# Reads all the PDB files from current directory and list them in an array
my $dir = getcwd;
my $dirName = basename($dir);
open(my $OUT, ">", "Redundant_".$dirName.".txt") or 
    die "Can not write file $!";
my @antibodyList = readDir ($dir);
@antibodyList = sort @antibodyList; 

# Variable Declaration
my (%redundantElements, %redundantClusters, %unique) = ();
my @identicals = ();
my ($ab1Light, $ab1Heavy, $ab2Light, $ab2Heavy);
my ($ab1Lstart, $ab1Lmid, $ab1Lend);
my ($ab1Hstart, $ab1Hmid, $ab1Hend);
my ($statusStartL, $statusEndL, $statusStartH, $statusEndH);
 
# Reads every element of the antibody list and compares with rest
for (my $ab1 = 0; $ab1 <= $#antibodyList ; $ab1++)
{
    # Checks if the antibody has already been compared and present in previous
    # cluster
    if (exists $redundantElements{$antibodyList[$ab1]})
    {
        next;
    }
    # Obtain Sequence (from PDB file) of each of antibody                     
    ($ab1Light, $ab1Heavy) = getAntibodySeq($antibodyList[$ab1]);
     # Split the Light and Heavy chains in 3 sections                          
#    next if ( ($ab1Light ) and (!$ab1Heavy) ) or
 #       ( (!$ab1Light ) and ($ab1Heavy) );
    
      

    ($ab1Lstart, $ab1Lmid, $ab1Lend) = splitString($ab1Light);    
    ($ab1Hstart, $ab1Hmid, $ab1Hend) = splitString($ab1Heavy);

    # Inner loop to compare the rest of element with the outer loop elemenet
    for (my $ab2 = ($ab1+1) ; $ab2 <= $#antibodyList; ++$ab2)
    {
	# Checks if the antibody has already been compared
	if (exists $redundantElements{$antibodyList[$ab2]})
	{
	    next;
	}
	# Obtain Sequence (from PDB file) of each of antibody 
	($ab2Light, $ab2Heavy) = getAntibodySeq($antibodyList[$ab2]);
	# Check if Middle section of ab1 is substring of ab2 
	# Then compare ends
	if ( (index($ab1Lmid, $ab2Light) != -1) or 
	     (index($ab2Light, $ab1Lmid)!=-1) ) # For light chain
	{
	    ($statusStartL, $statusEndL) = 
		isIdenticalEnds ($ab1Lmid, $ab2Light, $ab1Lstart, $ab1Lend); 
	}
	else
	{
	    next;
	}
	if ( (index($ab1Hmid, $ab2Heavy)!=-1) or 
	    (index($ab2Heavy, $ab1Hmid)!=-1) ) # For heavy chain
	{
	    ($statusStartH, $statusEndH) =
		isIdenticalEnds ($ab1Hmid, $ab2Heavy, $ab1Hstart, $ab1Hend);
        }
        else
	{
            next;
        }

        if ( ( $ab1Light) and ($ab1Heavy) )
            {
                # if all the 4 ends from both chains are identical then push them to
                # the identicals array (redundant)
                if (($statusStartL) && ($statusEndL) &&
                        ($statusStartH) && ($statusEndH))
                    {
                        push( @identicals, $antibodyList[$ab2] );
                        $redundantElements{$antibodyList[$ab2]} = 1; # Flag the redundants 
                    }
            }

        elsif ( ( $ab1Light) and (!$ab1Heavy)) {
            if ( ($statusStartL) && ($statusEndL) )
                {
                    push( @identicals, $antibodyList[$ab2] );
                    $redundantElements{$antibodyList[$ab2]} = 1;
                    # Flag the redundants  
                }
        }
        

        elsif ( ( $ab1Heavy) and (!$ab1Light)) {
            if ( ($statusStartH) && ($statusEndH) )
                {
                    push( @identicals, $antibodyList[$ab2] );
                    $redundantElements{$antibodyList[$ab2]} = 1;
                    # Flag the redundants  
                }
        }
        
        
    }
    


    

    

    if (@identicals)
    {                                         
	# hash with an array of redundant PDBs;     
	$redundantClusters{$antibodyList[$ab1]} = [@identicals];
	@identicals = ();                         
    }
    else 
    {                                                                        
	# hash with PDBs with no redundant PDB
	$unique{$antibodyList[$ab1]} = 1;     
	next;
    }

    
}

printClusters (\%redundantClusters, \%unique, $OUT);


=cut
=head1 NAME                                                                   


=head1 SYNOPSIS                                                               
use this program, type on command line,                                       
./getRedundantAntibodyClusters.pl

=head1 SYNOPSIS  
Perl Module for generating clusters for redundant antibodies                  

=head1 AUTHOR                                                                 
Saba Ferdous (ucbterd@acrm19)                                                 

=head1 COPYRIGHT AND LICENSE                                                  
Copyright (C) 2014 by Saba Ferdous                                            

This program is free software; you can redistribute it and/or modify          
it under the same terms as Perl itself, either Perl version 5.8.2 or,         
at your option, any later version of Perl 5 you may have available.           

=head1 BUGS                                                                   
                         
None reported... yet.

=cut




