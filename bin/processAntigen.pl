# This script runs on the directory of antibody-antigen complexes (processed 
# by get_complex_antibody program: antibody is numbered. PDB will not have any
# HETATM lines. These HETATMS are needed for antigens during epitope analysis)
# For that antigen is grabbed from original PDB where as antibody is just 
# copied to a new file.
# Now new PDB file called "pro^PDB.pdb" is formed with
# antigen from original PDB and antibody from the given PDB (processed). 

use strict; 
use warnings;
use general qw (readDirPDB); 
use epitope qw (getChainLabelAndType);
use pdb qw (get_pdb_path);
use SFPerlVars;

# Give directory address where antigens needs to process
my $dir = '/acrm/data/people/saba/data/dataNew/DataSJune2015/NR_Complex_Martin';
#my $dir = '/home/bsm/ucbterd/Desktop/NR_Antibody_Martin'; 
chdir $dir;

my @dirFiles = readDirPDB($dir);

foreach my $pdbFile (@dirFiles)
{
    chomp $pdbFile;
    # Get 1AFV from 1AFV_1.pdb
    my ($pdbCode, $ext) = split('_',  $pdbFile);
    
    # Get Antigen chain label
    my ($antigenChainLabel, $antigenChainType) =
	getChainLabelAndType($pdbFile);
   
    # Get original PDB path from local PDB
    my $pdbPath = get_pdb_path ($pdbCode);

   if ( ($antigenChainLabel eq "l") or ($antigenChainLabel eq "h"))
    {
        print "TEST: $pdbFile\n";
        
        #	$antigenChainLabel = "A"; # It was L, I replaced with A
        my ($antigenChainLabel, $antigenChainType) =
            getChainLabelAndType($pdbPath);
    }	

    # remove header first by program pdbatoms (returns only atom records)
    
    # Get Antigen chain and write on new file 
    my $pdbgetchain = $SFPerlVars::pdbgetchain;
    my $pdbatoms = $SFPerlVars::pdbatoms;
#    my @antigenPDB = `pdbatoms $pdbPath | $pdbgetchain $antigenChainLabel $pdbPath`;
    my @antigenPDB = `pdbatoms $pdbPath | $pdbgetchain $antigenChainLabel`; 
 #   print "@antigenPDB\n";
#exit; 
    
    if (!@antigenPDB)
    {
	$antigenChainLabel = "H";
	my @antigenPDB = `pdbatoms $pdbPath | $pdbgetchain $antigenChainLabel`;
    }

    # Open Output file in append mode
    my $outFile = "pro^$pdbFile";
    open (my $OUT, ">>$outFile") or
	die ("Can not open\n");
    
    # Write Antigen atom and HETATM records on new file 
    foreach my $line (@antigenPDB)
    {
	print {$OUT} $line if ($line !~ /HOH|PO4/); 
    }
    
    # Write Light chain on new file
    my @antibodyLightChain = `pdbatoms $pdbFile | $pdbgetchain L`;
    print {$OUT} @antibodyLightChain;
    
    # Write Heavy chain on new file
    my @antibodyHeavyChain = `pdbatoms $pdbFile | $pdbgetchain H`;
    print {$OUT} @antibodyHeavyChain;

    # Remove Old files 
    `rm $pdbFile`;
  #  last; 

}

# All the new files now have extension pro^:
# To remove that and restore original name, run the following bash command


# `for file in * ; do mv -v "$file" "${file#*^}"; done;`

