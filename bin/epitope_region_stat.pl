#!/usr/bin/perl -s 

# This Scripts reads all the antibody complexes from current/given directory
# and returns 2 files; one with the statistic of regions and odd bits and other 
# with the detailed information about epitopic sequence information (where regions
# and odd bits are separated by colon and multiple epitopic regions and odd bits
# are further separated by commas)
# Subroutines in module, Antigen_Contacts.pm 
# USAGE: ./epitope_region_stat.pl
# Author: Saba Ferdous

use strict;
#use warnings;
use Data::Dumper;
use Antigen_Contacts qw (antibody_antigen_contacts antibody_cont_residue      
antigen_cont_residue output_File_name get_hash_key get_fragments
get_regions_and_oddbits);
use general qw (readDirPDB);
#my ($gap, $contacting_residues);
#if ( defined ($::$gap) and defined ($::$contacting_residues) )
#{

print "Enter the allowed Gap\n";
my $gap = <STDIN>;
chomp $gap;

print "Enter the allowed contacting residues\n";
my $contacting_residues = <STDIN>;
chomp $contacting_residues;

# Reading direcoty files into an array     
#   my $dir = ".";                     
my $dir = "/acrm/data/people/saba/data/dataNew/DataMay2015/NR_Complex_Martin/ProteinAntigen";
chdir $dir; 
my @dir_files = readDirPDB($dir);
@dir_files = sort @dir_files;

# Open 2 files for writing, 1) For stats (Number of regions and odd bits) 
#2) For detail of residues
open(my $STAT, '>', "$dir/epitope-stat-G$gap-CR$contacting_residues") or die "File can not open";
open(my $EPITOPE_REGIONS, '>', "$dir/epitope_sequence-G$gap-CR$contacting_residues")
    or die "File Can not open";

print {$STAT} "Antibody:\tRegions:\tOdd_Bits\n";
print {$EPITOPE_REGIONS} "Antibody:Regions:Odd_Bits\n";

foreach my $pdb_file (@dir_files) 
{

    my ($antigen_chain_label, $antigen_resi, $antigen_chain_conts)
	= antigen_cont_residue($pdb_file);
    
    my @antigen_reseq = get_hash_key($antigen_chain_conts);
	
    my @fragments = get_fragments(\@antigen_reseq, $gap);
 
    my ($regions, $odds, $count_regions, $count_odd_bits)
	= get_regions_and_oddbits( $contacting_residues, @fragments);
    # To strip first comma - if present   
    $odds =~s/^,*//;
    print {$EPITOPE_REGIONS} "$pdb_file:$regions:$odds\n";
    print {$STAT} "$pdb_file:$count_regions:$count_odd_bits\n";   
	


}


mkdir -p "stats";
`mv ./epitope* ./stats`;




