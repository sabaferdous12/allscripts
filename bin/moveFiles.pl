# Move the files in list (input file) on a new directory
 
use strict;
use warnings;
use File::Copy;

my $file = $ARGV[0]; # List of files to be moved
my $destDir = $ARGV[1]; # Directory name to be moved in  

`mkdir -p $destDir`;

open(my $in_file, $file) or die "Can not open\n";

my @nr_files = <$in_file>;

print @nr_files;

#exit;
foreach my $nr_file (@nr_files){
    chomp $nr_file;
    
    if (-e $nr_file) {    
	move ($nr_file, $destDir);
    }

    else{
	print "$nr_file does not exist\n";
    }

}
