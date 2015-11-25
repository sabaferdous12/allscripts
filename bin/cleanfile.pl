# This short script removes spaces feom the output file of
# epitope statistics file
use strict; 
use warnings; 
my $infile = $ARGV[0]; 
open (my $IN, $infile);
open (my $OUT, ">$infile.Pro");
while (my $line = <$IN>)
{
   my ($pdb, $regions, $fragments) =
       split(/\s+/, $line);
   
   $pdb =~ s/^\s+//;
   $regions =~ s/^\s+//;
   $fragments=~ s/^\s+//;
   print {$OUT} $pdb,":", $regions,":", $fragments, "\n";
}
