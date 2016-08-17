use strict;
use warnings;
use Data::Dumper;


use general qw (readDirPDB);
my $inputFile = $ARGV[0];

my @PDBCodes = readFileDataInArray($inputFile);
use antibodyProcessing qw (
                              getChainTypeWithChainIDs
                              getPDBPath
                              getChains
                              readFileDataInArray
                              dirOperations
                      );
use antibodyAntigen qw (
                           getInterchainContacts
                           
                   );

my $ab = "L";

foreach my $pdb(@PDBCodes ) {
    chomp $pdb;
    
    print "Working on $pdb........\n";
    
    my $pdbPath = getPDBPath($pdb);
    my ($chainType_HRef, $chainIdChainTpye_HRef) =
        getChainTypeWithChainIDs ($pdbPath);
    my %chainType = %{$chainType_HRef};
    my @singleChainAb;
    
    
    if ( $ab eq "L") {
        @singleChainAb = @{$chainType{Light}};
    }
    elsif ( $ab eq "H") {
        @singleChainAb = @{$chainType{Heavy}};
    }
    
    my $pair;
    
    my %chainContacts = getInterchainContacts($pdbPath);
    print "INTER_CHAIN:", Dumper (\%chainContacts);
    
    my %dimers;
    
        foreach my $ab1 (@singleChainAb ) {
            
            foreach my $ab2 (@singleChainAb) {
                $pair = $ab1.$ab2;
                if ( exists $chainContacts{$pair}) {
   #                 if ( $chainContacts{$pair} > 80 ) {
                        $dimers{$pair} = $chainContacts{$pair};
 #                   }
                    
                }
            }
        }
    %dimers = reverse %{ { reverse %dimers } };
    
    print "$pdb....", Dumper (\%dimers);
    #last;
    
}

