#!/usr/bin/perl -s 
# This Scripts reads all the antibody complexes from current/given directory
# and returns 2 files; one with the statistic of regions and fragments and other 
# with the detailed information about epitopic sequence information (where regions
# and fragments are separated by colon and multiple epitopic regions and fragments
# are further separated by commas)
# Subroutines in module, Antigen_Contacts.pm 
# USAGE: ./epitope_region_stat.pl
# Author: Saba Ferdous

use strict;
#use warnings;
use Data::Dumper;
use Cwd;
use Antigen_Contacts qw
    (antibody_antigen_contacts
     antibody_cont_residue      
     antigen_cont_residue
     output_File_name
     get_hash_key
     getregionsAndfragments
     get_regions_and_oddbits
);
use general qw (readDirPDB);
print "Enter the allowed Gap\n";
my $gap = <STDIN>;
chomp $gap;

print "Enter the allowed contacting residues\n";
my $contacting_residues = <STDIN>;
chomp $contacting_residues;

# Reading direcoty files into an array     
my $dir = getcwd();
chdir $dir; 
my @dir_files = readDirPDB($dir);
@dir_files = sort @dir_files;

# Open 2 files for writing, 1) For stats (Number of regions and odd bits) 
#2) For detail of residues
open(my $STAT, '>', "$dir/epitope-stat-G$gap-CR$contacting_residues") or die "File can not open";
open(my $EPITOPE_REGIONS, '>', "$dir/epitope_sequence-G$gap-CR$contacting_residues")
    or die "File Can not open";
open (my $NONANTIGEN, '>', "$dir/epitope_nonAntigenPDBs") or die "Can not open file $!";

print {$STAT} "Antibody:Regions:Fragments\n";
print {$EPITOPE_REGIONS} "Antibody:Regions:Fragments\n";

foreach my $pdb_file (@dir_files) 
{
    print "TEST: $pdb_file\n";
    
    # if there are no contacts on any one of antibody chain (light or heavy)
    # then do not check antigen for epitopes - exclude non-antigens
    my ($light_cont_res_REFA, $heavy_cont_res_REFA, $light_chain_conts_REFH,
        $heavy_chain_conts_REFH) = antibody_cont_residue ($pdb_file);

    
    if ( (!@{$light_cont_res_REFA}) or (!@{$heavy_cont_res_REFA}) ) {
#        print {$EPITOPE_REGIONS} "$pdb_file:Not Antigen\n";
#        print {$STAT} "$pdb_file:Not Antigen\n";
        print {$NONANTIGEN} "$pdb_file\n";
        
        next;
    }
    
    else {
        my ($antigen_chain_label, $antigen_resi, $antigen_chain_conts, $antigen_Chains_HREF)
            = antigen_cont_residue($pdb_file);
        
        my @antigen_reseq = get_hash_key($antigen_chain_conts);
        	
        my @fragments = getregionsAndfragments(\@antigen_reseq, $gap, $antigen_Chains_HREF);
 
        my ($regions, $odds, $count_regions, $count_odd_bits)
            = get_regions_and_oddbits( $contacting_residues, @fragments);
        # To strip first comma - if present   
        $odds =~s/^,*//;
        print {$EPITOPE_REGIONS} "$pdb_file:$regions:$odds\n";
        print {$STAT} "$pdb_file:$count_regions:$count_odd_bits\n";   
	
    }
    
#    last; 
}


`mkdir -p "stats"`;
`mv ./epitope* ./stats`;




