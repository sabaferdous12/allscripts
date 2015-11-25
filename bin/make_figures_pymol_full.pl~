#!/usr/bin/perl
# This script reads epitope_sequence.txt file and makes figures to show regions
# and fragments in the 3D structure of an antigen by using pymol. It writes 
# .pml (pymol script) with all the commands and then sends that pymol script to
# program pymol. 

use strict;
#use warnings;
use Data::Dumper;
use SFPerlVars;
use epitope qw(getRegionRange getFragmentResidue writepymolScript );
use general qw (readDirPDB);

my @antigen_DNA;
my @antigen_protein;
my $aligned_pdb ;
my ($antigen_chain_label, $antigen_chian_type);
my (@rangeRegion, @residueFragment) = ();

my $dir = "/acrm/data/people/saba/data/dataNew/DataMay2015/NR_Complex_Martin";
chdir $dir; 
my $infile = "epitope_sequence-G3-CR3";
my $AntigenLength = "AntigenLength";

open (my $AG, '>', "./stats/$AntigenLength") or die
    "Can not open file\n";

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
	($antigen_chain_label, $antigen_chian_type) = 
	    split(" ", `$chaintype $pdb_file | head -1`);

        #### get Antigen length
        getAntigenChainLength($pdb_file, $antigen_chain_label, $AG);
                
        # To align the antibody structure around a center
	my $abalign = $SFPerlVars::abalign;
	$aligned_pdb = "aligned_".$pdb_file;
	`$abalign -k $pdb_file $aligned_pdb`;  
	
	@rangeRegion = getRegionRange($regions);
	@residueFragment = getFragmentResidue($fragments);
	
	    writepymolScript(\@rangeRegion, \@residueFragment, 
			     $antigen_chain_label, $aligned_pdb, $pdb_file);
	
    }
     

#    last;
} # While Loop ends here

# Directory Manipulation 

if (-d "Figures")
{
  if (-d "$infile")
  {
    `mv *.png ./Figures/$infile`;
  }
  else 
  {
   mkdir "./Figures/$infile";
   `mv *.png ./Figures/$infile`;
  }

}

else
{
mkdir "Figures";
if (-d "$infile")
  {
    `mv *.png ./Figures/$infile`;
  }
  else 
  {
   mkdir "./Figures/$infile";
   `mv *.png ./Figures/$infile`;
  }

}


print scalar @antigen_DNA, " DNA Antigens are:  @antigen_DNA\n";
print scalar @antigen_protein, " Protein Antigens are: @antigen_protein\n";
 


# ***************************Sub Routines ***********

sub getAntigenChainLength
    {
        my ($pdbFile, $antigenChain, $AG) = @_;
        my @chianData = split (" ", `pdbgetchain $antigenChain $pdbFile | pdbcount`);
        print {$AG} "$pdbFile: ", $chianData[3], "\n";
        
    }
    
