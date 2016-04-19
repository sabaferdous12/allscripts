use strict;
use warnings;

use general qw (readDirPDB);
my $dir = ".";

my @PDBfiles = readDirPDB($dir);
@PDBfiles = sort{ substr($a, 13) <=> substr($b,13) } @PDBfiles ;

print join ("\n", @PDBfiles );
#exit;


foreach my $pdb ( @PDBfiles) {
    chomp $pdb;

    print "$pdb\n";
    
    `pdbsecstr -s $pdb | awk '{""; print \$3}' >"$pdb.txt"`;
   # last;
#    `paste -d "" "$pdb.txt" >>ss. txt`;

    
}

#`for ( x in $(seq 1 1000); do 
 #   paste `echo -n "file$x "` > ss.txt

  #  done`


