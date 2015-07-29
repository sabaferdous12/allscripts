# Updated version
# This Script provides the clusters of redundant antibodies
# Version - 2
# It works by taking substring of protein (Middle part) being compared and then
# compares ends for 100% identity

use strict;
#use warnings;
use Data::Dumper;
use Cwd;
use File::Basename;
use redund qw(getAntibodySeq getFileData aa3to1 readDir splitString
getIndex compareStarts compareEnds checkSubstring);

# Reads all the PDB files from current directory and list them in an array
my $dir = getcwd;
my $dirName = basename($dir);
open(OUT, ">", "Redundant_".$dirName.".txt") or die "Can not write file $!";
my @antibodyList = readDir ($dir);

# Variable Declaration
my (%redundantElements, %redundantClusters, %unique) = ();
my @identicals = ();
my ($ab1Light, $ab1Heavy, $ab2Light, $ab2Heavy);
my ($ab1Lstart, $ab1Lmid, $ab1Lend);
my ($ab1Hstart, $ab1Hmid, $ab1Hend);
my ($ab2Lstart, $ab2Lend, $index_s);
my ($ab2Hstart, $ab2Hmid, $ab2Hend);
my ($statusStartL, $statusEndL, $statusStartH, $statusEndH);


# Reads every element of the antibody list and compares
for (my $ab1 = 0; $ab1 <= $#antibodyList ; $ab1++)
{
    # Checks if the antibody has already been compared
    if (exists $redundantElements{$antibodyList[$ab1]})
    {
        next;
    }
    # Inner loop to compare every element with the outer loop elemenet
    for (my $ab2 = ($ab1+1) ; $ab2 <= $#antibodyList; ++$ab2)
    {
	# Checks if the antibody has already been compared
	if (exists $redundantElements{$antibodyList[$ab2]})
	{
	    next;
	}

	# Obtain Sequence (from PDB file) of each of antibody 
	($ab1Light, $ab1Heavy) = getAntibodySeq($antibodyList[$ab1]);
	($ab2Light, $ab2Heavy) = getAntibodySeq($antibodyList[$ab2]);

	# Split the Light and Heavy chains in 3 sections
	($ab1Lstart, $ab1Lmid, $ab1Lend) = splitString($ab1Light);
	($ab1Hstart, $ab1Hmid, $ab1Hend) = splitString($ab1Heavy);

	# Check if Middle section of ab1 is substring of ab2 
	# Then compare ends
	if (index($ab1Lmid, $ab2Light) != -1)
	{
	    # Obtain start and ends of ab2
	    ($ab2Lstart, $ab2Lend, $index_s) = getIndex($ab2Light, $ab1Lmid);
	    # Compare the starts and ends and returns 0 or 1
	    $statusStartL = compareStarts($ab1Lstart, $ab2Lstart, $index_s);
   	    $statusEndL = compareEnds($ab1Lend, $ab2Lend);
	}

	elsif(index($ab2Light, $ab1Lmid)!=-1)
	{
	    ($ab2Lstart, $ab2Lend, $index_s) = getIndex($ab2Light, $ab1Lmid);
	    $statusStartL = compareStarts($ab1Lstart, $ab2Lstart, $index_s);
	    $statusEndL = compareEnds($ab1Lend, $ab2Lend);
	}

	else
	{
	    next;
	}
	
	if (index($ab2Hmid, $ab2Heavy)!=-1)
	{
            ($ab2Hstart, $ab2Hend, $index_s) = getIndex($ab2Heavy, $ab2Hmid);
            $statusStartH = compareStarts($ab1Hstart, $ab2Hstart, $index_s);
            $statusEndH = compareEnds($ab1Hend, $ab2Hend);
        }
	
	elsif(index($ab2Heavy, $ab2Hmid)!=-1)
	{
            ($ab2Hstart, $ab2Hend, $index_s) = getIndex($ab2Heavy, $ab2Hmid);
            $statusStartH = compareStarts($ab1Hstart, $ab2Hstart, $index_s);
            $statusEndH = compareEnds($ab1Hend, $ab2Hend);
        }
	
        else
	{
            next;
        }
	
	if (($statusStartL==1) && ($statusEndL==1) &&
	    ($statusStartH==1) && ($statusEndH==1))
	{
	    push( @identicals, $antibodyList[$ab2] );
	    $redundantElements{$antibodyList[$ab2]} = 1;
	}
	
    }

    if (@identicals)
    {                                         
	# An array is used as hash value;     
	$redundantClusters{$antibodyList[$ab1]} = [@identicals];
	@identicals = ();                         
    }
    
    else 
    {                                                                         
	$unique{$antibodyList[$ab1]} = 1;                                      
	next;
    }     
}



my (@modifiedCluster, @modifiedCluster2);

foreach my $key (sort {$a<=>$b} keys %redundantClusters)
{
    #Since hash value is an array - So it needs to dereferenced as below: 
    my @redundCluster = @{$redundantClusters{$key}};
    #The key and value are assigned to an array
    my @modifiedCluster = [$key, join(" ",@redundCluster)];
    my @modifiedCluster2 = split ("\n", "@{$modifiedCluster[0]}");
    
    # Each file is with.pdb extension, To avoid that in final output file, 
    # REGEX is used as below:

    foreach my $line (@modifiedCluster2)
    {
	$line =~ s/.pdb/,/g;
	print OUT "$line\n";
    }
}

foreach my $uniq (sort {$a<=>$b} keys %unique)
{
    $uniq =~ s/.pdb//g;
    print OUT "$uniq\n"; 
}


=head1 NAME                                                                   
get_redund.pl 

=head1 SYNOPSIS                                                               
use this program, type on command line,                                       
./get_redund.pl

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




