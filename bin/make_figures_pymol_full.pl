#!/usr/bin/perl
# This script reads epitope_sequence.txt file and makes figures to show regions
# and fragments in the 3D structure of an antigen by using pymol. It writes 
# .pml (pymol script) with all the commands and then sends that pymol script to
# program pymol. 

use strict;
#use warnings;
use Data::Dumper;
use SFPerlVars;
use epitope qw(getRegionRange getFragmentResidue);
use general qw (readDirPDB);
use antigenProcessing qw (getAntigenChains);

use Cwd;
use File::Basename;

my $pymol = "/usr/bin/pymol";

my @antigen_DNA;
my @antigen_protein;
my $aligned_pdb ;
my ($antigen_chain_label, $antigen_chian_type);
my (@rangeRegion, @residueFragment) = ();

my $dir = getcwd();
chdir $dir; 
my $infile = "epitope_sequence-G3-CR3";

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
    print "Processing ... $pdbFile\n";
    
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
	
	
#	my $chaintype = $SFPerlVars::chaintype;
        # To obtain the antigen chain label and chain type (N Protein)
#	($antigen_chain_label, $antigen_chian_type) = 
#	    split(" ", `$chaintype $pdb_file | tail -1`);

        my @antigenchains = getAntigenChains($pdb_file);
        
        # To align the antibody structure around a center
	my $abalign = $SFPerlVars::abalign;
	$aligned_pdb = "aligned_".$pdb_file;
	`$abalign -k $pdb_file $aligned_pdb`;  
	
	@rangeRegion = getRegionRange($regions);
	@residueFragment = getFragmentResidue($fragments);
	
	    writepymolScript(\@rangeRegion, \@residueFragment, 
			    \@antigenchains, $aligned_pdb, $pdb_file);
	
    }
     

#    last;
} # While Loop ends here

# Directory Manipulation 
`mkdir ComplexFigures`;
`mv *.png ComplexFigures`;


sub writepymolScript
    {
        my ($rangeRegionRef, $residueFragmentRef,
            $antigenChainLabel, $alignedPdb, $pdbFile) = @_;

        my $pdbFileName = basename($pdbFile, ".pdb");
        my @antigenChains = @{$antigenChainLabel};
                
        # Writing pymol script to provide as input to pymol program
        open (my $PYMOL, '>', "pymolscript.pml") or die
            "Can not open file";
        # opening perl script file
        
print {$PYMOL} <<__EOF;
load $alignedPdb
bg_color white
turn x, 90
turn y, 90
turn x, 90
turn y, 90
show cartoon
hide lines
color yellow, chain L
color pink, chain H
__EOF
        my $count = 1;
        foreach my $agChain (@antigenChains) {
            chomp $agChain;
            if ( ( $agChain eq "L") or ( $agChain eq "H") ) {
                $agChain = $count;
                print {$PYMOL} "color cyan, chain $agChain\n";
                $count++;
            }
            else {
                print {$PYMOL} "color cyan, chain $agChain\n";
            }
            
        }

       # Select the given epitope range and color that in red
        foreach  my $record (@{$rangeRegionRef})
            {
                my $chainLabel = substr ($record->[0], 0, 1);  # Obtaining first character
                my $start = substr ($record->[0], 1);  # Obtaining RESEQ - start
                my $end = substr ($record->[1], 1);  # Obtaining RESEQ - end

                my $count = 1;
                if ( ( $chainLabel eq "L") or ( $chainLabel eq "H") ) {
                    $chainLabel = $count;
                    print {$PYMOL} "select regions, resi $start-$end and chain $chainLabel\n";
                    print {$PYMOL} "color red, regions\n";
                    $count++;
                }
                else {
                    print {$PYMOL} "select regions, resi $start-$end and chain $chainLabel\n";
                    print {$PYMOL} "color red, regions\n";
                }

            }
        
        # Select and color the fragment residues
        foreach my $record (@{$residueFragmentRef})
            {
                my $chainLabel = substr ($record, 0, 1);
                 # Obtaining first character
                my $start = substr ($record, 1);

                my $count = 1;
                if ( ( $chainLabel eq "L") or ( $chainLabel eq "H") ) {
                    $chainLabel = $count;
                    print {$PYMOL} "select fragments, resi $start and chain $chainLabel\n";
                    print {$PYMOL} "color green, fragments\n";
                    $count++;
                }
                else {
                    print {$PYMOL} "select fragments, resi $start and chain $chainLabel\n";
                    print {$PYMOL} "color green, fragments\n";
                    
                }
                
            }

        # Taking the image
        print {$PYMOL} "ray 800,600\n";
        print {$PYMOL} "png $pdbFileName.png\n";
        print {$PYMOL} "quit\n";

        # Sending pymol script to program pymol
        `$pymol -c pymolscript.pml`;
        unlink "pymolscript.pml";
        unlink $alignedPdb;

}
