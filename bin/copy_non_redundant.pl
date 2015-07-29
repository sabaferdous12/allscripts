# Small script to take first element from file, which contains
# list of redundanat pdb file, in order to generate the list of
# non-redundant pdb files.

use strict;
use warnings;
use File::Copy;
#use Cwd; 

#my $dir = getcwd;
my $file = $ARGV[0];
my $out_dest = $ARGV[1];

open(my $in_file, $file) or die "Can not open\n";
my @nr_files;
while (<$in_file>){

my ($first_elem) = split(',', $_);
chomp $first_elem;
 
my $first_elem_pdb =  $first_elem.".pdb";
push (@nr_files, $first_elem_pdb);

}

#print join("\n", @nr_files);
#print scalar @nr_files;

foreach my $nr_file (@nr_files){
    chomp $nr_file;
    my $dest = "../$out_dest";
    if (-e $nr_file){
        copy ($nr_file, $dest);
    }

    else{
        print "$nr_file does not exists\n"; 
    }
}
