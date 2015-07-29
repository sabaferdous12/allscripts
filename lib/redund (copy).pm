package redund;

use Carp;

use strict;
#use warnings;

use Data::Dumper;

use Exporter qw (import);
our @EXPORT_OK = qw (get_seq_antibody read_dir get_file_data aa3to1 split_string get_index compare_strt compare_ends check_substring);

sub get_seq_antibody{

    my ($input_pro) = @_;
    
    my @protein = get_file_data($input_pro);
    
    my (@light,@heavy);
    foreach my $line (@protein){
	if ($line =~ /^ATOM*/){	
	    my @pro = split (' ', $line);
	    
	    if ($pro[4] eq 'L' and $pro[2] eq 'CA'){
		my $lg = aa3to1($pro[3]);
		push (@light, $lg);
	    }
	    
	    elsif($pro[4] eq 'H' and $pro[2] eq 'CA'){
		my $hv = aa3to1($pro[3]);
		push (@heavy, $hv);
	    
	    }
	    
	    
	}
	
    }
    my $seq_l = join('',@light);
    my $seq_h = join('',@heavy);
    
    
    return ($seq_l, $seq_h);
    

}


sub read_dir{

    my ($dir) = @_;
    my @list;
    opendir (DIR, "$dir") or die "Can't open directory: $!\n";
    while ( my $file = readdir(DIR)){
	next unless (-f "$dir/$file"); ## Reads only files
	
	next unless ($file =~ m/\.pdb$/); ## Reads files ending with pdb
	push(@list, $file);
	    }
	
	    closedir (DIR);
	    

	    return @list;
	}


sub get_file_data {
    my($filename) = @_;
  
    chomp $filename;
# Initialize variables                                                                                                
    my @filedata = ( );
    unless( open(FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }
    @filedata = <FILE_DATA>;
    close FILE_DATA;

    return @filedata;

}


sub aa3to1 {

    my($input) = @_;
    
    my %three2one = (
	'ALA' => 'A',
	'VAL' => 'V',
	'LEU' => 'L',
	'ILE' => 'I',
	'PRO' => 'P',
	'TRP' => 'W',
	'PHE' => 'F',
	'MET' => 'M',
	'GLY' => 'G',
	'SER' => 'S',
	'THR' => 'T',
	'TYR' => 'Y',
	'CYS' => 'C',
	'ASN' => 'N',
	'GLN' => 'Q',
	'LYS' => 'K',
	'ARG' => 'R',
	'HIS' => 'H',
	'ASP' => 'D',
	'GLU' => 'E',
	);

    # clean up the input
    $input =~ s/\n/ /g;
    
    my $seq = '';
    
    # This use of split separates on any contiguous whitespace
    my @code3 = split(' ', $input);
    
    foreach my $code (@code3) {
        # A little error checking
        if(not defined $three2one{$code}) {
            print "Code $code not defined\n";
            next;
        }
        $seq .= $three2one{$code};
    }
    return $seq;
}






sub split_string{

    my ($input_str) = @_;
    my ($strt_str, $mid_str, $end_str2, $temp_len, $temp_str);


    $strt_str = substr($input_str, 0, 5);
    $temp_str = substr($input_str, 5);
    $mid_str = substr($temp_str, 0, length($temp_str)-5);
#    $temp_len2 = length($input_str)-10;
 #   $mid_str = substr($input_str, 5, $temp_len2);
    $temp_len = length($strt_str) + length($mid_str);
    $end_str2 = substr($input_str, $temp_len);
    
    return ($strt_str, $mid_str, $end_str2);



}


sub get_index{

    my ($full_str, $sub_str) = @_;
    my ($strt_part, $end_part, $index_e, $index_s);

    $index_s = index($full_str, $sub_str);
    $strt_part = substr($full_str, 0, $index_s);
    $index_e = length($strt_part) + length($sub_str);
    $end_part = substr($full_str, $index_e);

    return ($strt_part, $end_part, $index_s);
}

sub compare_strt{

    my ($str1,$str2, $index) = @_;
    #my $temp = length($str1);
    my $count = 0;
    my $match = '';
    for (my $i = 4; $i>=0; $i--){
	if ($str1 && $str2){
	    if (substr($str1, $i, 1) eq substr($str2, $index-1, 1) ){
	    
		$match .= substr($str2, $index-1, 1);
		$count++;
		$index--;
	    }
	    else{
		next;
	    }
	}
       
    }
    my $rev_match = reverse($match);
   # print "$rev_match\n";
    my $status = check_substring($str1, $rev_match);

    return $status;
 
}

sub compare_ends{
    my ($str1, $str2) = @_;
    my $match = '';
    for (my $i = 0; $i<=4; $i++){
	if(!$str2){
	    if (substr($str1, $i, 1) eq substr($str2, $i, 1)){
		$match .= substr($str2, $i, 1);
	    }
	    else{
		next;
	    }
	}
	else{
	    next;}

 
    }

    #print "$match\n";
    my $status = check_substring($str1, $match);

    return $status;



}

sub check_substring{

    my ($str1, $str2) = @_;
    my $status;
    if (index($str1, $str2) != -1){
	$status = 1;
    }

    else {
	$status= 0;
    }
    return $status;
}


