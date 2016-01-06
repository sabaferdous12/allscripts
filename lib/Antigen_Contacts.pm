package Antigen_Contacts;

use 5.010001;
use strict;
#use warnings;
use SFPerlVars;
use antigenProcessing qw (getAntigenChains);

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(antibody_antigen_contacts get_antigen_chain_id 
antibody_cont_residue get_antibody_res_label get_antigen_res_label 
antigen_cont_residue output_File_name read_dir get_hash_key get_fragments
get_regions_and_oddbits);

our $VERSION = '0.01';


# Preloaded methods go here.

# Inputs: An antibody complex (pdb file)
# Outputs: Array containg contacts of antibody chains (LH) with antigen
#         e.g: Chain: L Res:  30A - Chain: A Res:  81  Contacts: 5
# Subroutine call: antibody_antigen_contacts($pdb_file); 
# Testing: my @chain_conts = antibody_antigen_contacts($pdb_file);
# Date: 02 April 2014
# Author: Saba
sub antibody_antigen_contacts{
    my ($pdb_file) = @_;
    my $chaincontacts = $SFPerlVars::chaincontacts;
    #my $antigen = get_antigen_chain_id($pdb_file);
    my (@chain_conts, @temp) ;
    
    my @antigenChains = getAntigenChains($pdb_file);
    foreach my $ag (@antigenChains ) {
        if ( $ag eq "L" or $ag eq "H") {
            $ag = "1";
        }
        
            print "TEST: $ag\n";
        
        chomp $ag;
        @chain_conts = `$chaincontacts -r 4.00 -x LH -y $ag $pdb_file`;
        splice @chain_conts, 0, 8;
        push (@temp, @chain_conts);
        
        #last;
        
    }
    
    
    print "TEST: @temp\n";
    exit;
    
    return @chain_conts;
}


# Inputs: An antibody complex (pdb file)                          
# Outputs: 2 Arrays of anonymous arrays, each with resSeq of light and heavy 
#          chains and number of contacts with  antigen. e.g: 30A 5
#          2 Hashes, each with ResSeq and total number of conatcts with antigen
#          Returning array and hash references
# Testing: scalar variables for keeping references
#          my ($light_conts, $heavy_conts, $light_chain_conts, 
#          $heavy_chain_conts) = &antibody_cont_residue($pdb_file); 
# Date: 02 April 2014                                                          
# Author: Saba
sub antibody_cont_residue{
    my ($pdb_file) = @_;
my    @chain_conts = antibody_antigen_contacts($pdb_file);
#    my (@chain_conts) = @_;
    my (@line_tokens, $chain_label, $res_num, @light_cont_res,
	@heavy_cont_res, $n);
    my $contacts = 0;
    my (%light_chain_conts, %heavy_chain_conts);
    foreach my $line (@chain_conts){
        chomp ($line);
	if ($line =~ /^Chain*/){
            @line_tokens = split /:/, $line;
            $chain_label = substr($line_tokens[1], 1, 1);
	    $res_num = get_antibody_res_label(@line_tokens);
	    $n = $line_tokens[5];
	    $n =~ s/^\s*(.*?)\s*$/$1/;

	    if ($chain_label eq 'L'){
	
		push(@light_cont_res, [$res_num, $n] );
		if(exists $light_chain_conts{$res_num}){
		    $contacts = $light_chain_conts{$res_num} + $n;
		    $light_chain_conts{$res_num} = $contacts;
		}
		else{
		    $contacts = $n;
		    $light_chain_conts{$res_num} = $contacts;
		}
		
            }elsif ($chain_label eq 'H'){
		push(@heavy_cont_res, [$res_num, $n] );
		if(exists $heavy_chain_conts{$res_num}){
		    $contacts = $heavy_chain_conts{$res_num} + $n;
		    $heavy_chain_conts{$res_num} = $contacts;
		}
		else{
		    $contacts = $n;
		    $heavy_chain_conts{$res_num} = $contacts;
		}
	    }

	}
    }

    return (\@light_cont_res, \@heavy_cont_res, \%light_chain_conts,
	    \%heavy_chain_conts);
}

# Inputs: An antibody complex (pdb file)                                      
# Outputs: Antigen chain label, an arrays of anonymous arrays, each with 
#          resSeq of light and heavy chains and number of contacts with 
#          antigen.
#          e.g: 30A 5                
#          A Hashes, each with ResSeq and total number of conatcts with antige 
#          Retuns array and hash references                                    
# Testing: scalar variables for keeping references                             
#          my ($ag_label, $antigen_conts, $antigen_contsss)
#          = &antigen_cont_residue($pdb_file);
# Date: 02 April 2014                                   
# Author: Saba 

sub antigen_cont_residue
{
    my ($pdb_file) = @_;
    my    @chain_conts = antibody_antigen_contacts($pdb_file);
  

    my (@line_tokens, $ag_chain_label, $res_num,  @antigen_conts, $n);
    my %antigen_chain_conts;
    my $contacts = 0;  
   foreach my $line (@chain_conts){
	chomp ($line);
        if ($line =~ /^Chain*/)
	{
            @line_tokens = split /:/, $line;
            $ag_chain_label = substr($line_tokens[3], 1, 1);
	    $n = $line_tokens[5];
	    $n =~ s/^\s*(.*?)\s*$/$1/;
	    if($ag_chain_label){
		$res_num = get_antigen_res_label(@line_tokens);		
		# array of anonymous arrays
		push(@antigen_conts, [$res_num, $n] );
	    
		if(exists $antigen_chain_conts{$res_num}){
                    $contacts = $antigen_chain_conts{$res_num} + $n;
                    $antigen_chain_conts{$res_num} = $contacts;
                }

                else{
                    $contacts = $n;
                    $antigen_chain_conts{$res_num} = $contacts;
                }

	    }
	}
    }
    return ($ag_chain_label, \@antigen_conts, \%antigen_chain_conts);

}


# Inputs: An array spliting a line into tokens
#        e.g: Chain: L Res:  30A - Chain: A Res:  81  Contacts: 5  
# Outputs: Returns reSeq of antibody. e.g: 30A
# Subroutine call: get_antibody_res_label(@line_tokens1);
# Testing: my reseq_antibody = get_antibody_res_label(@line_tokens1); 
# Date: 02 April 2014                                                         
# Author: Saba

sub get_antibody_res_label{
    my(@line_tokens1) = @_;
    my ($chain_res1, $reseq_antibody);
    $chain_res1 = $line_tokens1[2];
    $chain_res1 =~ s/^\s*(.*?)\s*$/$1/;
    $chain_res1 =~ m/(\d+[A-Z]?)/;
    $reseq_antibody = $1;

    return $reseq_antibody;

}

# Inputs: An array spliting a line into tokens                                 
#        e.g: Chain: L Res:  30A - Chain: A Res:  81  Contacts: 5   
# Outputs: Returns reSeq of antigen . e.g: 81                        
# Subroutine call: get_antigen _res_label(@line_tokens1                        
# Testing: my reseq_antigen  = get_antigen_res_label(@line_tokens1);         
# Date: 02 April 2014                                 
# Author: Saba 

sub get_antigen_res_label{
    my (@line_tokens) = @_;
    my ($chain_res, $reseq_antigen);
    $chain_res = $line_tokens[4];
    $chain_res =~ s/^\s*(.*?)\s*$/$1/; #Remves leading or trailing spaces
    $chain_res =~ m/(\d+[A-Z]?)/; # Retaining Insertion code for antibody 
                                  # Reseq e.g: 100A
    $reseq_antigen = $1;

    return $reseq_antigen;

}


# Inputs: An antibody complex (pdb file)
# Outputs: antigen label from list of chains in pdb file/complex 
#          e.g; A from 'ALH' 
# Subroutine call: get_antigen_chain_id($pdb_file);
# Testing: my $antigen_chain = get_antigen_chain_id($pdb_file); 
# Date: 02 April 2014                                                          
# Author: Saba 
sub get_antigen_chain_id {

    my ($pdb_file) = @_;
# using Andrew program got chian list like 'ALH'                               
    my $getchainlist = $SFPerlVars::getchainlist;
    my $chain_list = `$getchainlist $pdb_file`;
    my $antigen_chain = substr($chain_list, 1, 1); 
    $antigen_chain =~ s/^\s*(.*?)\s*$/$1/;
    
    return $antigen_chain;

}



# Inputs: An antibody complex (pdb file) name           
# Outputs: Output file name, e.g; for antibody 1AFV_1.pdb, output file would be
#          1AFV_1_contacts.xml
# Subroutine call: output_File_name($pdb_file)
# Testing: my $antigen_chain = get_antigen_chain_id($pdb_file);               
# Date: 02 April 2014
# Author: Saba  
sub output_File_name{
    my ($pdb_file) = @_;
    $pdb_file =~ s{\.[^.]+$}{}; # removes .pdb file extension
    my $xml_file = $pdb_file."_contacts.xml";
    return $xml_file;
}



# Inputs: path of Directory
# Outputs: Array with all files with pdb extension only
# Subroutine call: read_dir($dir)
# Testing: my $dir_files = read_dir($dir)
# Date: 02 April 2014      
# Author: Saba 
sub read_dir{
    my ($dir) = @_;
    my @list;
    opendir (DIR, "$dir") or die "Can't open directory: $!\n";
    while ( my $file = readdir(DIR)){
	next unless (-f "$dir/$file"); # Reads only files
	next unless ($file =~ m/\.pdb$/); # Reads files ending with pdb
	push(@list, $file);
    }
    closedir (DIR);
        
    return @list;
}


# Inputs: A hash reference
# Outputs: sorted keys of a hash in an Array
# Subroutine call: get_hash_key($hash_ref)
# Testing:my @hash_keys = get_hash_key($hash_ref)
# Date: 09 April 2014                                                               
# Author: Saba 
sub get_hash_key{
    my ($hash_ref) = @_;
    my @keys = keys % { $hash_ref };
    @keys = sort{$a<=>$b} @keys;
    return @keys;
}

# Inputs: A hash reference                                                          
# Outputs: values of a hash in an Array                                             
# Subroutine call: get_hash_val($hash_ref)                                         
# Testing:my @hash_values = get_hash_val($hash_ref)                                 
# Date: 09 April 2014                                                               
# Author: Saba
sub get_hash_val{
    my ($hash_ref) = @_;
    my @vals = values % { $hash_ref};

    return @vals;
}

# Inputs: An array reference                                                        
# Outputs: An Array of anonymous arrays                                             
# Subroutine call: get_fragments(\@antigen_reseq)
# Testing:my @fragments = get_fragments(\@antigen_reseq)
# Date: 09 April 2014                                                               
# Author: Saba
sub get_fragments{
    my ($antigen_reseq_ref, $gap) = @_;

    my $current_frag_ref;
    my @fragments = ( [] ); # Array of anonymous arrays
	$gap = $gap +1;
    foreach my $reseq (@$antigen_reseq_ref) {
        $current_frag_ref = $fragments[-1]; # Returns index of last array element
        my $length = scalar @{$current_frag_ref};
	# To check the cosective $reseq potentially interuppted by 1. 
	# Its equivalent to x = y+1 or x = y+2 -> x <= y + 2

        if ($length == 0 or $reseq <= @{$current_frag_ref}[-1]+$gap) { 
            push(@{$current_frag_ref}, $reseq);
        } else {
            push(@fragments, [$reseq]);        }
    }

return @fragments; 
}

# Inputs: An array of anonymous arrays                                               
# Outputs: 4 things; 2 strings, one string containing regions seperated by commmas,
#          one string containing odd-bits. 2 counters informing number of regions 
#          and odd bits 
# Subroutine call: get_regions_and_oddbits(@fragments)                               
# Testing:  my ($regions, $odds, $count_regions, $count_odd_bits)
#           = get_regions_and_oddbits(@fragments);
# Date: 09 April 2014                                                                
# Author: Saba 
sub get_regions_and_oddbits{
    my ($contacting_residues, @fragments) = @_;
    my $count_regions = 0;
    my $count_odd_bits = 0;
    my $regions = "";
    my $odds = "";

    foreach my $fragment (@fragments){
       my @fragment = @{$fragment};
       
       # To sort the numeric and alpha numeric values of array @fragment
       @fragment = sort { $a <=> $b || $a cmp $b } @fragment;

	 if (@$fragment>= $contacting_residues) 
	{

            $count_regions++; 
            if ($count_regions > 1) { # If not  first fragment, add a comma and 
		                      #concatenate as string
                $regions .= ",@fragment";
            } else {
		$regions = "@fragment"; # First fragment
            }
# Region is defined as a fragment of at least three consective residues
# potentially interupted by one residue 
	}
        else{
            # if fragment is of 2 consective numbers, treat them as one frag
            # while non-consective numbers would be treated as 2 frags
            if ( @$fragment == 2) {
                if ( $fragment[0] == ($fragment[1] - 1) ) {
                    $count_odd_bits++;
                }

                else {
		@fragment = join (",",@fragment);
                    $count_odd_bits = $count_odd_bits + 2;
                }
            }
            else {
                $count_odd_bits++;
            }

            
            if ($count_odd_bits > 1) {
		$odds .= ",@fragment";
            } else {
		$odds = "@fragment";
            }
         $odds =~ s/:,/:/; 
	}
       
    }
    
    return ($regions, $odds, $count_regions, $count_odd_bits);

}



1;
__END__

=head1 NAME

Antigen_Contacts - 

=head1 SYNOPSIS

  use Antigen_Contacts;

=head1 DESCRIPTION
 This Perl module contains the several subroutine for calculation of contacts 
 between antigen and antibody and finding the epitopic regions, odd bits in contact with antibody

=head2 EXPORT

None by default.

=head1 AUTHOR

Saba Ferdous, E<lt>ucbterd@localdomainE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Saba Ferdous

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


