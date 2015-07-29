use strict;
use warnings;

my $dirPath1 = $ARGV[0];
my $dirPath2 = $ARGV[1];

open (my $FILE1, ">Combined.txt") or die "Can not open file";
open (my $FILE2, ">Merged.txt") or die "Can not open file";

my @Combined = read_dir($dirPath1);
my @Merged = read_dir($dirPath2);
print {$FILE1} join ("\n", @Combined);
print {$FILE2} join ("\n", @Merged);

#`comm -13 <(sort Combined.txt) <(sort Merged.txt) (>difference.txt)`;
#`grep -Fxvf ./Combined.txt ./Merged.txt >difference.txt`;



########################


sub read_dir{
    my ($dir) = @_;
    my @list;
    opendir (DIR, "$dir") or die "Can't open directory: $!\n";
    while ( my $file = readdir(DIR)){
	next unless (-f "$dir/$file"); # Reads only files
	next unless ($file =~ m/\.pdb$/); # Reads files ending with pdb
	my ($a, $b) = split('_', $file);
	push(@list, $a);
    }
    closedir (DIR);
        
    return @list;
}
