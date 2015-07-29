#!/usr/bin/perl -s 
# This small script reads xml file from command line and
# output the list of PDB codes on terminal which can be
# directed to a file

use strict; 
use lib("~/scripts/lib");
use general qw (getPdbCodeListFromXmlFile);

main(); 

sub main ()
{ 
UsageDie() if defined($::h);

my $SACsFile = $ARGV[0];
my @PDBsList = getPdbCodeListFromXmlFile ($SACsFile);
print join ("\n", @PDBsList), "\n";
}

sub UsageDie
{
    print <<__EOF;
    getSACsPDBList V1.1 (c) 2014, UCL, Saba Ferdous

	Usage: getSACsPDBList.pl <Input XML File>
__EOF
exit 0;
}
