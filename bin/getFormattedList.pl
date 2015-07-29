use strict;
use warnings;

my ($DIFF, $COMBINED, $OUT);
my (@free_antibody, @complex);
open ($DIFF, '<', "difference.txt") or die "Can not open file $DIFF\n";
open ($COMBINED, '<', "./Redundant_files/Redundant_Combined_Martin.txt")
    or die "Can not open file $COMBINED\n";

open ( $OUT, '>', "Free_Complex.txt")
  or die "Can not open file\n";

print $OUT "Free Antibody:Complex\n";

my @combined = <$COMBINED>;

while (my $line = <$DIFF>)
{
    
    chomp $line;
    my ($records) = grep (m/$line/, @combined);
    my @redundants = split (/,\s+/, $records);
    
    foreach my $elem (@redundants)
    {
	chomp $elem;
	my $elem_pdb = $elem.".pdb";
	if (-e "./Complex_Martin/$elem_pdb")
	{
	    push (@complex, $elem);
	}
	else
	{
	    push (@free_antibody, $elem);
	}
	next if ( (!@complex ) or (!@free_antibody) );

    }
    print $OUT join(',', @free_antibody), ":", join(',', @complex), "\n";
    @free_antibody = ();
    @complex = ();
}

