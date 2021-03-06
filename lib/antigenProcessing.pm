package antigenProcessing;
use strict;
use warnings;

use Exporter qw (import);
our @EXPORT_OK = qw (
                        getAntigenChains
                );

sub getAntigenChains
{
    my ($pdbFile) = @_;
    my (@antigenChainsTemp, @antigenChains, @antigens) ;

    # Grep all antigen chain in PDB file
    @antigenChainsTemp = `grep "CHAIN A" $pdbFile`;
        
    # For multiple antigens in PDB file
    if ( (scalar @antigenChainsTemp) > 1  )
    {
        foreach my $ag ( @antigenChainsTemp )
        {
            chomp $ag;
            @antigens = split(m/\s+/, $ag);
            push (@antigenChains, $antigens[5]);
        }
    }
    # For Single antigen
    else {
        @antigens = map {split(" ", $_)} @antigenChainsTemp;
        push (@antigenChains, $antigens[5]);
    }
    
    return @antigenChains;

}
