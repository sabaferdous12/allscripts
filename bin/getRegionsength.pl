use strict; 
use warnings; 
use Data::Dumper; 

use lib("~/allscript/lib");
use epitope qw(getRegionRange); 

open (my $IN, "./stats/epitope_sequence-G3-CR3") or die $!; 
open (my $OUTR, ">RegionLengthTEST") or die "Can not open file";
open (my $OUTF, ">FragmentLengthTEST")or die "Can not open file";

my @file = <$IN>;
my @rangeRegion; 

foreach my $line (@file)
{
    next  if ($line =~ /^Antibody/); 
    my ($pdb_file, $regions, $fragments) = split (':', $line);
    chomp $fragments;

    # To obtain the range on each epitopic region i.e, start and end                                         
    @rangeRegion = getRegionRange($regions);

    # To find the length of each region
    my $rangeLen; 
	print {$OUTR} $pdb_file, ": ";	
    foreach my $range (@rangeRegion) 
        {
            $rangeLen = ($$range[1]-$$range[0])+1;
            print {$OUTR} "$rangeLen "; 
        }
	print {$OUTR} "\n";

    # To find the length of each fragment
    if ( $fragments ) {
        my @frags = split (",", $fragments);
        my (@fragSize, @fragLen);
print {$OUTF} $pdb_file, ": ";
        for ( @frags ) {
            chomp $_;
            @fragLen = split (" ", $_);
            print {$OUTF} scalar @fragLen, " ";
        }
	print {$OUTF} "\n";
    }
    else {
	print {$OUTF} $pdb_file, ": ";	
	print {$OUTF} 0, "\n"; 
        next;
    }
}
