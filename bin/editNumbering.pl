#!/usr/bin/perl -s
# This script reads a numbered ab-ag complex and renumbers antibody
# by chosed numbering scheme.
# usage: ./editNumbering <numbering scheme> <Dir Name>
#        For example: ./editNumbering -c . 
use strict;
use warnings;
use SFPerlVars;
use general qw (readDirPDB);

my ($nsch, $numbering);
my $kabatnum;
$kabatnum = $SFPerlVars::kabatnum;
 
if ( $::k) {
    $nsch = '-k';
    $numbering = "Kabat";    
}

if ( $::c) {
    $nsch = '-c';
    $numbering = "Chothia";    
}

if ( $::a) {
    $nsch = '-a';
    $numbering = "Martin";    
}
# Reading directory into an array
my $inDir = $ARGV[0];
my @dirFile = readDirPDB($inDir);

# Loop throug each file in the dircetory
foreach my $inFile (@dirFile )
{
    print "Renumbering: $inFile\n";
    my $abFile = "tempAB.pdb";
    my $headers = "headers_Ag.pdb";
    
    open (my $IN, '<', $inFile) or die "Error - Cannot open input file $!\n ";
    open (my $AB, '>', $abFile) or die "Error - Cannot open file for writing\n";
    open (my $HEAD_AG, '>', $headers) or die "Error - Cannot open file for writing\n";
    
    # Write antibody chains in this file
    while ( my $line = <$IN>) {
        if ( ( ( $line =~ /\s+L\s+/) and ($line =~/^ATOM/) )
                 or ( ( $line =~ /\s+H\s+/) and ($line =~/^ATOM/) ) ) {
            print {$AB} $line;
        }
        else {
            # write headers and antigen information in this file
            print {$HEAD_AG} $line;        
        }
        
    }
    close $HEAD_AG;
    close $AB;
    close $IN;
    
    # Apply numbering
    `$kabatnum $nsch $abFile numbered.pdb`;
    
    open ($HEAD_AG, '<', $headers) or die ":(";
    open (my $NUM, '<', "numbered.pdb");
    open (my $FINAL, '>', "final.pdb");
    
    # Read header_Antigen file for only headers
    # Write them on new file final.pdb
    while ( my $header =  <$HEAD_AG> )
    {
        if ( $header =~ /^REMARK/)
        {
            # Replace numbering name
            if ( $header =~ /NUMBERING/) {
                ($header) =~ s/$header/REMARK 950 NUMBERING  $numbering\n/;
                print {$FINAL} $header;
            }
            else {
                print {$FINAL} $header;
            }
        }
    }
    close $HEAD_AG;

    # Read numbered.pdb and write on final.pdb
    while ( my $num = <$NUM>)
    {
        next if ($num =~ /^MASTER|^END/);
        print {$FINAL} $num;
    }
    close $NUM;
    
    open ($HEAD_AG, '<', $headers) or die ":(";
    # Read header_antigen file for only antigen
    while (  my $header =  <$HEAD_AG>)
    {
        if ( ( $header =~ /^ATOM/) or ( $header =~ /^HETATM/) ) {
            print {$FINAL} $header;
        }
    }

    close $HEAD_AG;
    close $FINAL;
    
    unlink $inFile;
    unlink $abFile;
    unlink $headers;
    unlink "numbered.pdb";

    `mv final.pdb $inFile`;

} # End iteration for each file

