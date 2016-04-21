#!/usr/bin/perl
# This script reads epitope_sequence.txt file and makes figures to show regions
# and fragments in the 3D structure of an antigen by using pymol. It writes 
# .pml (pymol script) with all the commands and then sends that pymol script to
# program pymol. 
# It also gives a text file with length of antigen chain
use strict;
#use warnings;
use Data::Dumper;
use SFPerlVars;
use epitope qw(getRegionRange getFragmentResidue writepymolScript );
use general qw (readDirPDB);
use antigenProcessing qw (getAntigenChains);

use Cwd;
 
my @antigen_DNA;
my @antigen_protein;
my $aligned_pdb ;
my ($antigen_chain_label, $antigen_chian_type);
my (@rangeRegion, @residueFragment) = ();

my $dir = getcwd();

chdir $dir; 
my $infile = "epitope_sequence-G3-CR3";

open (my $AGPEP, '>', "./stats/peptideAntigenLength") or die
    "Can Not Open File...\n";

open (my $AGPRO, '>', "./stats/proteinAntigenLength") or die
    "Can Not Open File\n";

open (my $AGPP, '>', "./stats/peptideAntigenPDBs") or die
    "Can Not Open File\n";

open (my $AGPR, '>', "./stats/proteinAntigenPDBs") or die
    "Can Not Open File\n";


my @dirFiles = readDirPDB ($dir); 
my $EPITOPE; 

# Open file containg Epitope sequence
open ($EPITOPE, '<', "./stats/$infile") or die
    "Can not open file";

my @epitopeSequenceFile = <$EPITOPE>;


# Reading every PDB file in the directory
foreach my $pdbFile (@dirFiles)
{
    chomp $pdbFile; 
    print "Processing... $pdbFile\n";
    my $epitopeInfo; 
    # Looking for Epitope Sequence in the file
    foreach my $epitope (@epitopeSequenceFile)
    {
	chomp $epitope; 
	
	if ($epitope =~ /$pdbFile/)
	{
	    $epitopeInfo = $epitope;  
	}
    }
    
    # Obtaining PDB file name, regions and fragments in variable    
    my ($pdb_file, $regions, $fragments) = split (':', $epitopeInfo);
    
    if ($_ =~ /^Antibody/) # ignoring first line
    {
	next;
    }
    
    else
    {
	open (my $COMPLEX, $pdb_file) or die
	    "Can not open $pdb_file"; # open PDB file
	
	
	my $chaintype = $SFPerlVars::chaintype;
        # To obtain the antigen chain label and chain type (N Protein)
#	($antigen_chain_label, $antigen_chian_type) = 
#	    split(" ", `$chaintype $pdb_file | tail -1`);
        my @antigenChains = getAntigenChains($pdb_file);
                
        #### get Antigen length
        #getAntigenChainLength($pdb_file, $antigen_chain_label, $AGPEP, $AGPRO, $AGPP, $AGPR);
                
        # To align the antibody structure around a center
	my $abalign = $SFPerlVars::abalign;
	$aligned_pdb = "aligned_".$pdb_file;
	`$abalign -k $pdb_file $aligned_pdb`;  
	
	@rangeRegion = getRegionRange($regions);
	@residueFragment = getFragmentResidue($fragments);

        
	    writepymolScript(\@rangeRegion, \@residueFragment, 
			     \@antigenChains, $aligned_pdb, $pdb_file);
        
    }
     

#    last;
} # While Loop ends here

# Directory Manipulation 
`mkdir -p Figure`;
`mv *.png ./Figure`;

print scalar @antigen_DNA, " DNA Antigens are:  @antigen_DNA\n";
print scalar @antigen_protein, " Protein Antigens are: @antigen_protein\n";
 


# ***************************Sub Routines ***********

sub getAntigenChainLength
    {
        my ($pdbFile, $antigenChain, $AGPEP, $AGPRO) = @_;
        my @chainData = split (" ", `pdbgetchain $antigenChain $pdbFile | pdbcount`);

        if ( $chainData[3] < 30 ) {
            print {$AGPEP} "$pdbFile: ", $chainData[3], "\n";
            print {$AGPP} "$pdbFile\n";
        }
        else {
            print {$AGPRO} "$pdbFile: ", $chainData[3], "\n";
            print {$AGPR} "$pdbFile\n";
        }
        
    }
    
