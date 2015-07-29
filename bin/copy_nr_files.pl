# Copy the files in list (input file) on a new directory
 
use strict;
use warnings;
use File::Copy;

my $file = "/acrm/data/people/saba/data/dataNew/DataMay2015/NR_Complex_Martin/stats/peptideAntigenPDB";

open(my $in_file, $file) or die "Can not open\n";

my @nr_files = <$in_file>;

#print @nr_files;

#exit;
foreach my $nr_file (@nr_files){
    chomp $nr_file;
    my $dest = "/acrm/data/people/saba/data/dataNew/DataMay2015/NR_Complex_Martin/PeptideAntigen";
    if (-e $nr_file){    
	move ($nr_file, $dest);
    }

    else{
	print $nr_file;
    }

}
