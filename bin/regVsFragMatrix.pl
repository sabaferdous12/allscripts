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
open (my $IN, $infile) or die "Can not opn file\n";
my @epitopeStats =  <$IN>;


my $csvFile = 'regionFragmentRes.csv';
my $csv = Text::CSV->new;
open (my $OUT, '>:encoding(utf8)', $csvFile) or
    die "Could not open $csvFile: $!";

my $count = 0;
my @regionFragmentCount2d;
my @regionFragmentCount = ();  

# Loop throughs the regions 
# Loop goes upto the maximum regions

my @region = qw (
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
45
46
47
50
51
80
);

my @frags = qw (0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
);


for (my $i = 0; $i < 46; $i++)
#foreach my $f (@frags)
{
    @regionFragmentCount = (); # Resetting the array
    # Loop throughs the fragments     
    # Loop goes upto the maximum fragments
    for (my $j = 0; $j < 17; $j++)
        {
	# Reads every line to count the number of regions against the number 
	# of fragments and forms a array of arrays
	foreach my $epitope (@epitopeStats)
	{    
	    chomp $epitope;
	      # Split the data line on colon
            my ($pdbID, $regions, $fragments) =
		    split(/:/, $epitope);
#		print "SSSS: :::: $pdbID: $regions: $fragments\n";		

		    if (($regions == $region[$i]) and ($fragments == $frags[$j])) 
                    {
			$count++; 
                    }
            
	}

	push (@regionFragmentCount, $count);
        
	$count = 0;
    }
    # Constructing array of arrays
    push (@regionFragmentCount2d, [@regionFragmentCount]);
  #  last;
}
# For printing the array of arrays
foreach my $regionFragmentCount (@regionFragmentCount2d)
{
    my @regionFragmentCount = @{$regionFragmentCount}; 
 #   print join (" ", @regionFragmentCount), "\n"; 
    $csv->print ($OUT, $regionFragmentCount);
    print $OUT "\n"; 
}

