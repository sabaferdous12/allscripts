#!/usr/bin/perl -w
# movePDBs.pl --- struct;
# Author: Saba Ferdous <ucbterd@martin-cd01.biochem.ucl.ac.uk>
# Created: 20 Apr 2016
# Version: 0.01

# Fir given list of PDB codes, moves all the PDB file containing the PDB code
# e.g: For 1CIC it will move 1CIS_1 and 1CIS_2 ..... 
use warnings;
use strict;

use File::Copy;
my $file = $ARGV[0]; # List of files to be moved
my $destDir = $ARGV[1]; # Directory name to be moved in

`mkdir -p $destDir`;

open(my $in_file, $file) or die "Can not open\n";
my $dir = ".";


my @inFiles = <$in_file>;


foreach my $infile( @inFiles) {
    chomp $infile;

    opendir ( my $DIR, $dir ) or
        die "Can not open directory $dir";
    
    foreach  my $fl ( grep m/$infile/, readdir ( $DIR ) )
    {
        my $from = $dir."/".$fl;
        move $from, $destDir;
        print "MOVED: $infile\n";
    }
}

__END__

=head1 NAME

movePDBs.pl - Describe the usage of script briefly

=head1 SYNOPSIS

movePDBs.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for movePDBs.pl, 

=head1 AUTHOR

Saba Ferdous, E<lt>ucbterd@martin-cd01.biochem.ucl.ac.ukE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2016 by Saba Ferdous

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
