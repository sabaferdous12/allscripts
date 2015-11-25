# This script reads epitope statistics file,
# constructs a 2d matrix (array of arrays) by counting 
# the number of fragments againt number of regions and stores results 
# in a csv file as a matrix  
use strict;
use warnings; 
use Data::Dumper; 
use Text::CSV; 
#use Data::Dump qw (dump);
use Array::Transpose; 

my $infile = $ARGV[0];
open (my $IN, $infile);
my @epitopeStats =  <$IN>;

my $csvFile = $infile.".csv";
my $csv = Text::CSV->new;
open (my $OUT, '>:encoding(utf8)', $csvFile) or
    die "Could not open $csvFile: $!";

my $count = 0;
my @regionFragmentCount2d;
my @regionFragmentCount = ();  

# Loop throughs the regions 
# Loop goes upto the maximum regions
my @regions = qw (0 1 2 3 4 5 6 7 8 9);
my @frags = qw (0 1 2 3 4 5 6 7 8 9 10 11 12 13 14);
for (my $i = 0; $i < 10; $i++)
{
    @regionFragmentCount = (); # Resetting the array
    
    # Loop throughs the fragments     
    # Loop goes upto the maximum fragments
    for (my $j = 0; $j < 15; $j++)
    {
	# Reads every line to count the number of regions against the number 
	# of fragments and forms a array of arrays
	foreach my $epitope (@epitopeStats)
	{	    
	  # Split the data line on colon
	    my ($pdbID, $regions, $fragments) =
	    split(":", $epitope);

	    if (($regions == $regions[$i]) and ($fragments == $frags[$j])) 
	    {
		$count++; 
	    }
	}
	push (@regionFragmentCount, $count);
	$count = 0;
    }
    # Constructing array of arrays
    push (@regionFragmentCount2d, [@regionFragmentCount]);
}

# For printing the array of arrays
foreach my $regionFragmentCount (@regionFragmentCount2d)
{
    my @regionFragmentCount = @{$regionFragmentCount}; 
 #   print join (" ", @regionFragmentCount), "\n"; 
    $csv->print ($OUT, $regionFragmentCount);
    print $OUT "\n"; 
}

