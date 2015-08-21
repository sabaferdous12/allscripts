#! /usr/bin/perl
# This short script will check the type of chain in given PDB (Antibody-
# Antigen Complex) and will identify the Antigen as protein or DNA or RNA.
# Will create a directory for non-protein antigens and will move non-protein
# antigen in that directory


use strict; 
use warnings; 
use general qw(readDirPDB); 
use SFPerlVars;
use Cwd;

# Read directory of Antibody-Antigen Complexes
# my $dir = "/acrm/data/people/saba/data/dataNew/DataMar2015/NR_Complex_Martin";
my $dir = getcwd();

chdir $dir; 
my @dirFiles = readDirPDB($dir);

# Open files for writing list of nonprotein and protein antigens
open(my $PRO, '>', "ProteinAntigen.list") or 
    die ("Can not open file $!"); 

open(my $NONPRO, '>', "nonProteinAntigen.list") or
    die ("Can not open file $!");
#print "@dir_files\n"; 

`mkdir -p nonProteinAntigen`;
 
my (@nonProteinAntigen, @proteinAntigen);

# reading each PDB file in the directory 
foreach my $pdbFile (@dirFiles)
{    
    open (my $COMPLEX, $pdbFile) or die
	"Can not open $pdbFile"; # open PDB file
    
    my $chainType = $SFPerlVars::chaintype;
        # To obtain the antigen chain label and chain type (N Protein)
    my  ($antigenChainLabel, $antigenChianType) = 
	split(" ", `$chainType $pdbFile | head -1`);

    if ( ($antigenChianType eq 'DNA') or 
	 (($antigenChianType eq 'RNA') ) )
        # To keep track of DNa antigens
    {
	push (@nonProteinAntigen, $pdbFile);
	`mv $pdbFile nonProteinAntigen`; 	
    }
    else
    {
	push (@proteinAntigen, $pdbFile);
	next; 
    }

}

print {$PRO} join ("\n", @proteinAntigen);
print {$NONPRO} join ("\n", @nonProteinAntigen);

`mv ProteinAntigen.list nonProteinAntigen`;
`mv nonProteinAntigen.list nonProteinAntigen`;
