use strict;
use warnings;
use Data::Dumper;
use Array::Transpose;

my @array;

my %resHash;
my $count = 1;
my @matrix;

#print"Residue, Coil, Helix, Strand";
push ( @matrix, ["Residue", "Coil", "Helix", "Strand"]);

while (my $line=<>) {
    chomp($line);
    #print length ($line);
#  print "Residue: Coil, B-bridge, Bend, Turn, A-Helix, 5-Helix, 3-Helix\n";  
    if ( $line =~ m/^"\S*"/) #and (length ($line) > 1000) ){
       {
#           print length ($line);
    
#        print $line, "\n";
        $resHash{$count} = $line;
        
        my ( $coilCount, $alphaCount, $strandCount ) = getSScount ($line);
        push (@matrix, [$count, $coilCount, $alphaCount, $strandCount]);
        
        #print "$count: $coilCount, $strandCount, $alphaCount\n";
        
        $count++;
        
 #   last;    
       }
    
    #last;
    
}

#print Dumper (@matrix);

my $array=transpose(\@matrix);
foreach my $l ( @{$array} )
    {
        print join (",", @{$l} ) , "\n";
    }


sub getSScount
{
    my ($SS) = @_;
    my ($coilCount, $strandCount, $alphaCount ) =
            (0, 0, 0);
    my @ss = split ("", $SS);
    foreach my $ss (@ss )
    {
        if ( ( $ss eq "~" ) or ( $ss eq "-") or ( $ss eq "S") or ($ss eq "s")
         or ( $ss eq "T" ) or ( $ss eq "t") ) 
        {
            $coilCount++;
        }
        elsif ( ( $ss eq "B" ) or ( $ss eq "b") or ( $ss eq "E") or ( $ss eq "e") )
        {
            $strandCount++;
        }

        elsif ( ( $ss eq "H") or ($ss eq "h") or ($ss eq "G") or ($ss eq "g")
            or ( $ss eq "I" ) or ( $ss eq "i") )
        {
            $alphaCount++;
        }
 
    }
    my ($cper, $aper, $sper);
    $cper = (($coilCount*10)/1000000) * 100;
    $aper = (($alphaCount*10)/1000000) *100;
    $sper = (($strandCount*10)/1000000) *100;

    $cper = sprintf "%.2f", $cper;
    $aper = sprintf "%.2f", $aper;
    $sper = sprintf "%.2f", $sper;
        
#    print "Number Coil: $coilCount\n";
    return $cper, $aper, $sper;
    
    
}

#print Dumper (\%resHash);
