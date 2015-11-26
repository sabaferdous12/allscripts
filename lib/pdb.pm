package pdb;
#***********************************************************************
#
#   Perl Module:   pdb.pm
#   File:          pdb.pm 
#
#   Version:       V1.2
#   Date:          22.04.14
#   Fucntion:      The perl module is a library for program processAntibodyPDBs
#                  and contains most of the subroutines for data process and 
#                  uses external programs from ~/perl5/lib/perl5/SFPerlVars.pm
#   Modules:       SFPerlVars, general
#
#   Copyright:  (c) UCL, Saba Ferdous, 2014                                   
#   Author:     Miss Saba Ferdous                                             
#   Address:    Institute of Structural and Molecular Biology                 
#               Division of Biosciences                                       
#               University College                                            
#               Gower Street                                                  
#               London                                                        
#               WC1E 6BT                                                      
#   EMail:      saba@bioinf.org.uk  
#
#***********************************************************************
use Carp;
use IO::CaptureOutput qw ( capture capture_exec qxx qxy );
use strict;
use warnings;
use Data::Dumper;
use File::Copy;
use SFPerlVars;
use general qw (readDirPDB); 
use Exporter qw (import);
#use pdbWrap qw (movePDBs);

our @EXPORT_OK = qw (
    get_pdb_path
    get_file_data
    get_pdbcode_from_xml  
    check_chain 
    split_pdb_to_chains 
    get_pdbchains_contacts 
    largest_value_hash 
    pair_heavy_light_antigen 
    get_hash_key 
    get_hash_val 
    antibody_assembly 
    antibody_number 
    extract_CDR 
    assemble_CDR_antigen 
    antigen_CDR_conts 
    get_complex 
    get_antibody_antigen_complex 
    make_antibody_complex
    hasHapten
    processHapten
    getchainIDs
    getChainTypeWithChainIDs
    getSingleChainAntibodyAntigenComplex
    getResolInfo
    getsingleChainAntibody
    checkHashforIdenticalvalues
                   );

#use pdbWrap qw (movePDBs);

# ************ get_pdbcode_from_xml **************
# Description: Reads pdb codes from xml file (Obtained for SACS)
#              into an array and returns an array
# Inputs: XML file from SACS
# Outputs:An array with PDB codes
# Subroutine call/Testing: my @pdb_codes_XML = get_pdbcode_list_xml($xml_file);
# Date: 26 June 2014                                 
# Author: Saba

sub get_pdbcode_from_xml
{
    my ( $input_file ) = @_;

    open ( my $IN, "<$input_file" ) or die "Can not open $input_file\n";
    my @pdb_lines = 
	grep /pdb/, <$IN>; # Read file in array / grep line with word pdb 
    my @pdb_codes;

    foreach my $line ( @pdb_lines )
    {
	my ( $start, $code) = split ( "\"", $line );
	push ( @pdb_codes, $code );
    }

    close $IN;
    return @pdb_codes;
}

# ************* get_pdb_path *****************
# Description: For given pdb it return local pdb file path 
# Inputs: PDB code
# Outputs: Complete file path for local pdb e.g: /acrm/data/pdb/pdb3u7a.ent
# Subroutine call/Testing: $file_path = get_pdb_path($pdb_id)
# Date: 26 June 2014                                         
# Author: Saba 

sub get_pdb_path
{
    my ( $pdb_id ) = @_;
    chomp $pdb_id;
    my $file_path;

    if ( $pdb_id=~m/[0-9A-Z]{4}/ )
    {
	my $pdbprep = $SFPerlVars::pdbprep;
	my $pdbext = $SFPerlVars::pdbext;
	
	$file_path = $pdbprep . lc ( $pdb_id ) . $pdbext;
    }
    else
    {
	croak "This is not a valid PDB Code\n";
    }
    
    return $file_path;
   
}

# ************* get_file_data *****************                                
# Description: Reads a text file data into an array line by line
# Inputs: A text file with PDB code on each line
# Outputs: An array with PDB codes
# Subroutine call/Testing: @file_data = get_file_data ($txt_file)
# Date: 26 June 2014       
# Author: Saba
                                                                 
sub get_file_data 
{
    my ( $filename ) = @_;
    chomp $filename;
    my @filedata = ( );
    my $FILE_DATA;

    unless ( open ( $FILE_DATA, $filename ) )
    {
	croak "Cannot open file \"$filename\"\n\n";
	
    }
    @filedata = <$FILE_DATA>;
 
    close $FILE_DATA;

    return @filedata;

}


# ************* check_chain *****************                                
# Description: Checks chain type either Light. Heavy and Antigen 
# Inputs: PDB file / Here PDB file on a particular location
# Outputs: Three scalars; each for different chain type 
#          e.g: Light, Heavy and Antigen
# Subroutine call/Testing: my ($lg, $hv, $ag) = check_chain($file_path);
# Date: 26 June 2014                                                    
# Author: Saba  

sub check_chain
{
    my ( $pdb_path )=@_;
    my $countpdb = $SFPerlVars::countpdb;
    # Program pdbcounts number of chains, atoms and residues
    my $pdb_count = `$countpdb $pdb_path`;

    my $idabchain = $SFPerlVars::idabchain;
    my @chain_info = `$idabchain $pdb_path`;#Chain 1, A: Antigen    
    # Program idabchain identifies chain type and label

    my @chains = split ( ' ', $pdb_count );
    my $chain_no = $chains[1];
    my @chain_type;
    
    foreach my $line ( @chain_info )
        {
            my @chainInfo = split ( ' ', $line );
            push( @chain_type, $chainInfo[3] );
        }

    # Variable Declaration
    my ( $antignNum, $heavyNum, $lightNum ) = (0);
    my ( $antigen, $heavy, $light ) = ("");
    foreach my $rec ( @chain_type )
    {
	if ( grep $rec eq 'Antigen', @chain_type )
	{
	    $antignNum++;
	    $antigen = 'Antigen';
	}
	elsif ( grep $rec eq 'Heavy', @chain_type )
	{
	    $heavyNum++;
	    $heavy = 'Heavy';
	}
	elsif ( grep $rec eq 'Light', @chain_type )
	{
	    $lightNum++;
	    $light = 'Light';
	}

    }
    if ( ( !$light) and (!$heavy) and (!$antigen) )
        {
            return 
    }
    return ($light, $heavy, $antigen);

}

# ************* split_pdb_to_chains *****************
# Description: Splits all chains of pdb into separate files
# Inputs: PDB file / Here PDB file on a particular location                   
# Outputs: Returns number of chains - a PDB file has been splitted to
# Subroutine call/Testing: $chainno =  split_pdb_to_chains ($file_path);
# Date: 26 June 2014       
# Author: Saba

sub split_pdb_to_chains
{
    my $dir = ".";
    my ( $pdbfile ) = @_;
    open ( my $PDBFILE, $pdbfile) or die "Could not open file \"$pdbfile\"" ;
    # This program splits chains in a single PDB file into multiple PDBs
    my $pdbSplitChains = $SFPerlVars::pdbsplitchains;
    my $pdbatoms = $SFPerlVars::pdbatoms; 
    my $pdbgetchain = $SFPerlVars::pdbgetchain; 
    # Reading input file from stdin 
    #`$pdbSplitChains < $pdbfile` or cat $pdbfile | $pdbSplitChains;   
    # To strip header and keeping only ATOMs 
    `$pdbatoms $pdbfile | $pdbSplitChains`;

    my @chainPDBs = readDirPDB($dir); # Reads all PDB files in directory
    foreach my $chain (@chainPDBs)
    {
	my ($chainID, $ext) = split (/\./, $chain);
       # To discard HETATM by selecting each chain and keeping only ATOM lines
        `$pdbgetchain -a $chainID $chain >$chainID`;
        `mv $chainID $chainID.pdb`; # Rename A as A.pdb

    }
}


# ************* get_pdbchains_contacts *****************                      
# Description: Calculates the inter-chain contacts and return them into a hash
# Inputs: PDB file / Here PDB file on a particular location               
# Outputs: Returns a hash with chain pair as key and number of contacts as val
# Subroutine call/Testing: %hash = get_pdbchains_contacts($file_path);
# Date: 26 June 2014                                                          
# Author: Saba

sub get_pdbchains_contacts
{
    my ( $pdb_file ) = @_;
    # The program chaincontacts calculates inter chain contacts 
    my $chaincontacts = "$SFPerlVars::chaincontacts -r 4.00";
    my @chain_conts = `$chaincontacts $pdb_file`;
    
    # To chop the header of program output
    splice @chain_conts, 0, 8;
    my %chainTable=();
    my $ncontacts=0;

    foreach my $line ( @chain_conts )
    {
	chomp ($line);
	##Chain: A Res:  79  - Chain: H Res: 105  Contacts: 13
	if ($line =~ /^Chain*/)
	{
	    my @cont= split /:/, $line;
	    my $n = $cont[5];
	    chomp($n);
	    my $ch1 = $cont[1];
	    my $ch2 = $cont[3];
	    chomp($ch1);
	    chomp($ch2);
	    $ch1 =~ s/Res//;
	    $ch2 =~ s/Res//;
	    $ch1 =~ s/^\s*(.*?)\s*$/$1/; # Removing leading and trailing space
	    $ch2 =~ s/^\s*(.*?)\s*$/$1/;
	    # Create Chain's pair e.g: LH
	    my $chkey=$ch1.$ch2;
	    
	    # Counting number of contacts between 2 chains and storing them 
	    # in a hash like LH 150 
	    if ( exists $chainTable{$chkey} )
	    {
		$ncontacts = $chainTable{$chkey} + $n;
		$chainTable{$chkey}=$ncontacts;
	    }
	    else
	    {
		$ncontacts=$n;
		$chainTable{$chkey}=$ncontacts;
	    }
	}
    } 
    
    return %chainTable;   
    
}

# ************* largest_value_hash *****************                      
# Description: Finds the largest value in a hash 
# Inputs: Takes hash as input
# Outputs: returns key with largest value
# Subroutine call/Testing: my ($key, $val) = largest_value_hash(\%hash);
# Date: 26 June 2014                                                      
# Author: Saba                                                          

sub largest_value_hash 
{
    my $hash_ref   = shift;
    my ( $key, @keys ) = keys   %$hash_ref;
    my ( $big, @vals ) = values %$hash_ref;
    
    for ( 0 .. $#keys ) 
    {
        if ( $vals[$_] > $big ) 
	{
            $big = $vals[$_];
            $key = $keys[$_];
        }
    }
    return $key, $big;
}

# ************* pair_heavy_light_antigen *****************                     
# Description: This sub-routines does multiple jobs:
#              1. returns light and heavy pairs with inter-chain contacts 
#                 in a hash     
#              2. returns antigen inter-chain contacts in a hash
#              3. returns antigen labels 
# Inputs: PDB file / path to PDB file 
# Outputs: 2 hashes - one with antibody (light and heavy chains) with contacts
#                     one with antigen-antigen contacts
#          1 array - containg antigen chain labels
# Subroutine call/Testing: my ($heavy_light, $antigen_pair) = 
#                              &pair_heavy_light_antigen($file_path);
# Date: 26 June 2014                                                       
# Author: Saba

sub pair_heavy_light_antigen
{
    my ( $pdb_path )=@_;
    my ( @heavy, @light, @antigen, %heavy_light_pair, %ag_hash );
          
    my $idabchain = $SFPerlVars::idabchain;
    my @chain_info = `$idabchain $pdb_path`;

    my %chainsHash = getChainTypeWithChainIDs (@chain_info);
    @antigen = @{$chainsHash{Antigen}};
    @light =  @{$chainsHash{Light}};
    @heavy =  @{$chainsHash{Heavy}};
   
    
# Some of antigens have chain labels of H and L (Same as antibody chains). 
# These need to be renamed. The reason is to avoid the problem that arises
#  when grouping light and heavy LH with antigen to compute number of contacts
# in CDRs of antibody with antigen. 
# e.g: In cases with antigen labelled as H will give wrong results. LH->H
# Therefore, this bit of code is renaming the antigen by converting the 
# chain label to lower case and renaming the PDB file as %h or %l
# and updating the antigen labels array
  
    my $renumpdb = $SFPerlVars::renumpdb;
    my $index=0;

    foreach my $antigen_chain( @antigen )
    {
        
        chomp $antigen_chain;
        
	if( ( $antigen_chain eq "L" ) or ( $antigen_chain eq "H" ) )
	{
	    my $H_antigen = $antigen_chain.".pdb";
	    my $newChainLabel = lc ($antigen_chain);
    `$renumpdb -c $newChainLabel -n -d $H_antigen "%"$newChainLabel.pdb`;
	    my $filename = "%".$newChainLabel; 
	    rename ( $H_antigen, "Backup_$antigen_chain.pdb" );
           
            # Replacing H or L by % chain label
            
 	    splice (@antigen, $index, 1, $filename); # Replacing H with %h 
	#    push ( @antigen, $filename ); # Replacing H or L by % chain label
	    
	}

	$index++;
	
    }

    # using contact hash (hash with inter-chain contact information)
    my %contacts_hash = get_pdbchains_contacts(@_);
    my %new_hash = ();
 
    # find light and heavy pair with maximum contacts
    foreach my $lg ( @light )
    {
	foreach my $hv ( @heavy )
	{
	    my $pair = $lg.$hv;
	    if ( exists $contacts_hash{$pair} )
	    {
		$new_hash{$pair} = $contacts_hash{$pair};
	    }
	    elsif ( !exists $contacts_hash{$pair} )
	    {
		$new_hash{$pair} = 0;
	    }
	}
	# Select key with largest value 
	my ( $key, $val ) = largest_value_hash( \%new_hash );
	$heavy_light_pair{$key}= $val;
	%new_hash = (); # Empty hash
    }
    # find antigen-antigen chain contact pairs
    foreach my $ag1( @antigen )
    {
	foreach my $ag2( @antigen )
	{
	    my $ag_pair = $ag1.$ag2;
	    if ( exists $contacts_hash{$ag_pair} )
	    {
		$ag_hash{$ag_pair} = $contacts_hash{$ag_pair};
	    }
	    
	}
    }
    
    %ag_hash = reverse %{ {reverse %ag_hash} }; # removing duplicate values
    return ( \%heavy_light_pair, \%ag_hash, \@antigen, \%chainsHash);
}

# ************* get_hash_key *****************                        
# Description: Finds all keys of a hash and return into an array 
# Inputs: Takes hash as input                                              
# Outputs: Array with hash keys - here antibody chain labels e.g: LH MK
# Subroutine call/Testing: my @keys = get_hash_key($heavy_light); 
# Date: 26 June 2014                                                        
# Author: Saba

sub get_hash_key
{
    my ( $hash_ref ) = @_;
    my @keys = keys % { $hash_ref };
    
    return @keys;
}

# ************* get_hash_val *****************                            
# Description: Finds all values of a hash and return into an array             
# Inputs: Takes hash as input                                                
# Outputs: Array with hash values - here antibody chain contacts e.g: 198 186  
# Subroutine call/Testing: my @vals = get_hash_val($heavy_light);
# Date: 26 June 2014                                                      
# Author: Saba 

sub get_hash_val
{
    my ( $hash_ref ) = @_;
    my @vals = values % { $hash_ref};

    return @vals;
}

# ************* file_assembly *****************                              
# Description: It assembles two files together for a given chain labels pair. 
#              In our case it assembles L and H chain of an antibody in 1 PDB
# Inputs: An array contaning chain label pairs
# Outputs: It returns number of antibodies assembled and writes separate file 
#          for each
# Subroutine call/Testing: my $count = file_assembly (@keys);
# Date: 26 June 2014                                                          
# Author: Saba                                                                 

sub antibody_assembly
{
    my ( @chain_pair_label ) = @_;
    my $dir = '.';
    my $count = 0;
    my ($L_READ, $H_READ, $OUT);
    foreach my $pair ( @chain_pair_label )
    {
	my @chain = ();
	my ( $L, $H ) = split ( '', $pair );
	
	my $l_pdb = $L.".pdb";
	my $h_pdb = $H.".pdb";
	open ( $L_READ, '<', "$dir/$l_pdb" ) or
	    die "Can not open file $l_pdb\n";
	open ( $H_READ, '<', "$dir/$h_pdb" ) or
	    die "Can not open file $h_pdb\n";
	open ( $OUT, '>', "$dir/$pair.pdb" ) or
	    die "Can not write file $pair.pdb\n";
	
	while ( !eof ($L_READ ) )
	{
	    my $light = <$L_READ>;
	    #next if ($light =~ /^HETATM/);
	    print $OUT $light;
	    
	}
	
	while ( !eof ($H_READ ) )
	{
	    my $heavy = <$H_READ>;
	    #next if ($heavy =~ /^HETATM/);
	    print $OUT $heavy;
	    
	}
	
	$count++;
    }
    close $L_READ;
    close $H_READ;
    close $OUT;
        
    return $count;
}


# ************* antibody_number *****************                             
# Description: Numbers the antibody according to the user defined numbering 
#              scheme 
# Inputs: A reference to array containing antibody (LH) chain labels and
#         numbering scheme
# Outputs: It returns number of antibodies numbered and writes separate file   
#          for each numbered antibody  
# Subroutine call/Testing: my $count1 = antibody_number (\@keys, $nsch);
# Date: 26 June 2014      
# Author: Saba

sub antibody_number
{
    my ( $antibody_pairs, $nsch ) = @_;
    my ( $out, $error, $numbered_pdb );
    
    foreach my $pair( @$antibody_pairs )
    {
	my $antibody = $pair.".pdb";
	my $kabatnum = $SFPerlVars::kabatnum;
	$numbered_pdb = $pair."_num.pdb";

# Capturing error from Stderr
# In case of error (kabat program failure), output file ($numbered_pdb)
# would be empty 
	($out, $error) = qxx ( "$kabatnum $nsch $antibody $numbered_pdb" );
		
	if ( -z $numbered_pdb )
	{ # Checks if file is empty
            croak $error;
        }
    }

}

# ************* extract_CDR *****************                                
# Description:Extracts CDR regions from antibody according to kabat definition 
# Inputs: A reference to array containing antibody (LH) chain labels
# Outputs: returns count of processed antibodies and 2 files containing 
#          structure information of CRDs regions for light and heavy chains
# Subroutine call/Testing: my $count2 = extract_CDR(\@keys);
# Date: 26 June 2014                                                        
# Author: Saba

sub extract_CDR
{
    my ( $antibody_pairs, $ab ) = @_;
    my ( $out, $error, $success, $CDR );

    foreach my $pair ( @$antibody_pairs )
    {
	my $antibody = $pair."_num.pdb";
	my $getpdb= $SFPerlVars::getpdb;
      	my $CDR_pdb = $pair."_CDR.pdb";
        my @numbering;
        
        if ( $ab eq "LightHeavy") {
            @numbering =  (["L24", "L34"], ["L50", "L56"], ["L89", "L97"],
                              ["H31", "H35"], ["H50", "H65"], ["H95", "H102"] );
        }
        elsif ( $ab eq "Light") {
            @numbering =  (["L24", "L34"], ["L50", "L56"], ["L89", "L97"]);
        }

        elsif ( $ab eq "Heavy") {
            @numbering =  (["H31", "H35"], ["H50", "H65"], ["H95", "H102"] );
        }

        my $all_error = '';
	# to avoid append and clearing the file contents if already present
	if (-e $CDR_pdb)
	{	   
	    open ($CDR, '>', $CDR_pdb) or die "can not open $CDR_pdb";

	}	    

	foreach my $cdr (@numbering) 
	{
	    my ( $start, $end ) = @{$cdr};
	    ($out, $error) = 
		qxx ( "$getpdb $start $end $antibody >>$CDR_pdb" );
	    if ( $error )
	    {
		$all_error .= $error;
	    }
	    
	}
	    
	if ( $all_error )
	{
	    croak $all_error;
	}

    }
#    close ($CDR); 
}


# ************* assemble_CDR_antigen *****************                     
# Description: Assembles CDR regions with antigen and writes them into a file 
#              in order to calculate the number of contacts of each antibody 
#              (CRDs) with antigen 
# Inputs: 2 array references containing antibody (LH) chain labels and
#         antigen chain labels    
# Outputs:  Returns a hash of hash containg antibody contacts with 
#           correnponding antigen (computed based on number of contacts) 
#           between CDRs and antigen 
# Subroutine call/Testing: my %antigen_conts = 
#                          assemble_CDR_antigen(\@keys, $antigen_labels);
# Date: 26 June 2014                                                          
# Author: Saba 

sub assemble_CDR_antigen
{
    my ( $ab_pair, $antigen_chains ) = @_;
    my $dir = '.';
    my ( $AG_FILE, $CDR_FILE, $AG_CDR_FILE, $cdr_ag_conts );
    my %return_hash = ();
    my $nonAgCount;
    
    foreach my $antibdy ( @$ab_pair )
    {
	my $cdrs = $antibdy."_CDR.pdb";
	$return_hash{$antibdy} = {};

	foreach my $ch ( @$antigen_chains )
	{
	    my $pdb = $ch.".pdb";
	    next if ( $pdb eq '0.pdb');# to deal exception in 3J30 (split prob)
	    
	    open ( $AG_FILE, '<', "$dir/$pdb" ) or
		die "Could no open file $pdb\n";
	    open ( $CDR_FILE, '<', "$dir/$cdrs" ) or
		die "Could not open file $cdrs\n";
	    open ( $AG_CDR_FILE, '>',"$ch"."_"."$cdrs" ) or 
		die "Can not write file $cdrs\n";
	    
	    while ( !eof ( $AG_FILE ) )
	    {
		my $antigen_pdb = <$AG_FILE>;
		print $AG_CDR_FILE $antigen_pdb;
	    }
	    
	    while ( !eof ($CDR_FILE ) )
	    {
		my $cdrs_pdb = <$CDR_FILE>;
		print $AG_CDR_FILE $cdrs_pdb;
	    }
	    
	    my $ag_cdr_filename = "$ch"."_"."$cdrs";
# To calculate the contacts of antigen with CDRs of antibody
            my $cdr_ag_cont;
            
	     $cdr_ag_cont = antigen_CDR_conts ( $ag_cdr_filename, $ch );
	    $return_hash{$antibdy}->{$ch} = $cdr_ag_cont;
	    
	}
	
    }
    
    close ($AG_FILE);
    close ($CDR_FILE);
    close ($AG_CDR_FILE);  
    return  %return_hash; 
    
}

# ************* antigen_CDR_cont *****************                          
# Description: Computes number of contacts between CDRs and antigen
# Inputs: PDB file and antigen chain label 
# Outputs: A scalar containing number of contacts an antigen is making with 
#          antibody
# Subroutine call/Testing: my $no_contacts
#                               = antigen_CDR_conts($file_path, $antigen);
# Date: 26 June 2014      
# Author: Saba 

sub antigen_CDR_conts
{

    my ($pdb_file, $antigen_chain ) = @_;
    my $Lchaincontacts =
	"$SFPerlVars::chaincontacts -r 4.00 -x L -y $antigen_chain";
    my @Lcdr_conts = `$Lchaincontacts $pdb_file`;

  my $Hchaincontacts =
      "$SFPerlVars::chaincontacts -r 4.00 -x H -y $antigen_chain";
    my @Hcdr_conts = `$Hchaincontacts $pdb_file`;
    my $n;            
    splice @Lcdr_conts, 0, 8;
    splice @Hcdr_conts, 0, 8;
        
    # if any of antibody chain (L or H) doesn't have contacts with CDRs
    # Then do not compute contacts with antigen
    if ( (@Lcdr_conts) and (@Hcdr_conts) )
        {
            my @cdr_conts = (@Lcdr_conts, @Hcdr_conts);
            $n = countCDRContacts (@cdr_conts);
        }
    elsif ( (@Lcdr_conts) or (@Hcdr_conts) ) {
        my @cdr_conts = (@Lcdr_conts, @Hcdr_conts);
        $n = countCDRContacts (@cdr_conts);
            
        if ( $n >= 50 ) {
            $n = $n;
        }
        else {
            $n = 0;
        }
    }
    
    
    else {
        $n =0;
    }
    
    return $n;
    
}

sub countCDRContacts
    {
        my (@cdr_conts) = @_;
        my $n = 0;            
        foreach my $line ( @cdr_conts )
            {
                chomp ($line);
                
                if ( $line =~ /^Chain*/ )
                    {
                        my @conts = split /:/, $line;
                        $n += $conts[5];
                        chomp($n);
                        
                    }
            }
        return $n;
    }
    
# ************* antigen_CDR_cont *****************                            
# Description: Computes hash of antibody and antigen complex on basis of 
#              highest number of contacts
# Inputs: A hash with contact information of CDRs with each of antigen
# Outputs: A hash with antibody chain labels as key and associated antigen 
#          chain label as value
# Subroutine call/Testing: my %complex_hash = get_complex(%antigen_conts); 
# Date: 26 June 2014                                                         
# Author: Saba 

sub get_complex
{
    my ( %cdr_ag_contact_hash ) = @_;
    my %complex_hash;
       
    foreach my $ab_pair (keys %cdr_ag_contact_hash)
    {

            my %antigen_chains = %{ $cdr_ag_contact_hash{$ab_pair} };
	
            # Get highest value from hash and corrrespoding key 
            my $Largest_v =  ((sort {$a <=> $b} values %antigen_chains)[-1]);
            my $Largest_k = ((sort { $antigen_chains{$b} <=> $antigen_chains{$a} } keys %antigen_chains)[0]);
            
            # Get second highest value from hash and corrrespoding key 
            my $secLargest_v =  ((sort {$a <=> $b} values %antigen_chains)[-2]);
            my $secLargest_k = ((sort { $antigen_chains{$b} <=> $antigen_chains{$a} } keys %antigen_chains)[1]);
            
            # To check if second highest value is non-zero
            if ( $secLargest_v > 0 ) {
                $complex_hash{$ab_pair} = [$Largest_k, $secLargest_k];
            }
            
            elsif ( ( $Largest_v > 0 ) and ($secLargest_v == 0 ) ) 
                {
                    $complex_hash{$ab_pair} = [$Largest_k];
                }
            else
                {
                    $complex_hash{$ab_pair} ="NULL"; 
                }

    }
    

    return %complex_hash;
}

# ************* checkHashforIdenticalvalues *****************                
# Description: This subroutine checks if all the values of hash are "NULL"
#              If yes then returns 1 otherwise zero
# Inputs: A hash
# Outputs: Flag of 0 or 1
# Subroutine call/Testing: my $nonAgFlag =  checkHashforIdenticalvalues
#                                            (%complex_hash);
# Date: 01 Oct 2015
# Author: Saba  
sub checkHashforIdenticalvalues
    {
        my (%antibodyNonAg) = @_;
        my $hashSize = keys %antibodyNonAg;
        my $countNonAg = 0;
        my $nonAgFlag = 0;

        foreach my $ab (keys %antibodyNonAg ) {
            if ( $antibodyNonAg{$ab} eq "NULL") {
                $countNonAg++;
            }
        }
        
        if ( $countNonAg == $hashSize) {
            $nonAgFlag = 1;
        }
        else 
            {
                $nonAgFlag = 0;
            }
            return $nonAgFlag;

        }
# ************* get_antibody_antigen_complex *****************                
# Description: Assembles the antibody with corresponding antigen based on 
#              information in complex_hash 
# Inputs: PDB code of query PDB and complex hash
# Outputs: returns number of files processed and Writes complex files 
# Subroutine call/Testing: my $count3 = get_antibody_antigen_complex ($pdb_id,
#                                         %complex_hash);
# Date: 26 June 2014
# Author: Saba  

sub get_antibody_antigen_complex
{
    my ( $pdb_id, $destABAG, $destAB, $chainsHashRef, $numbering,
         $pdb_path, %complex_hash) = @_;
    my $dir = '.';
    my $count = 1;
    my $biAntigen = 0;
    my %chainsHash = %{$chainsHashRef};
    my %mapedChains;
    
    my @antibodychains = grep {$complex_hash{$_} eq "NULL" } keys %complex_hash;

   $count =  make_antibody_complex(\@antibodychains, \$pdb_id, "", $count,
                                   $numbering, $pdb_path, $chainsHashRef);
      
    movePDBs ($dir, $destAB, $pdb_id);
    my @antigen;
    
    foreach my $ab_pair( keys %complex_hash )
    {
######
        my $antibody = $ab_pair."_num.pdb"; 
	 my $antigenRef;
        
	next if ( $complex_hash{ $ab_pair } eq  "NULL" );

        if ( $complex_hash{ $ab_pair } ne  "NULL" )
	{
	    $antigenRef = $complex_hash{$ab_pair};
            @antigen = @{$antigenRef};
	} 
          
#####
        my ($L, $H) = split ("", $ab_pair);

#        my $Ag = $complex_hash{$ab_pair};
        
        my @LightChains = @{$chainsHash{"Light"}};
        my @HeavyChains = @{$chainsHash{"Heavy"}};
        my @AntigenChains = @{$chainsHash{"Antigen"}};        
        if ( grep (/$L/, @LightChains) ){
            $mapedChains{L} = $L;
        }
        elsif ( grep (/$L/, @HeavyChains)) {
            $mapedChains{H} = $L;
        }
        if ( grep (/$H/, @LightChains)) {
            $mapedChains{L} = $H;
        }
        elsif ( grep (/$H/, @HeavyChains) ) {
            $mapedChains{H} = $H;
        }
       # if ( grep (/$Ag/, @AntigenChains)) {
        #    $mapedChains{A} = $Ag;
        #}


        open ( my $AB_FILE, '<', "$dir/$antibody" ) or
            die "Could not open file $antibody";
        
        open ( my $AG_AB_FILE, '>>',"$pdb_id"."_"."$count".".pdb" ) or
            die "Can not write complex";

                
        my @Ag;
        
        foreach my $agChain( @antigen) {
            if ( grep (/$agChain/, @antigen)) {
                push (@Ag, $agChain);
            }
        }
        $mapedChains{A} = [@Ag];   
        printHeader($AG_AB_FILE, $numbering, $pdb_path, %mapedChains);    

        foreach my $agChain(@antigen) {
            
        $agChain = $agChain.".pdb";
           
        
######
	open ( my $AG_FILE , '<', "$dir/$agChain" ) or
	    die "Could no open file $agChain";
            
            while ( !eof ( $AG_FILE ))
                {
                    my $file1 = <$AG_FILE>;
                    # To ignore the Footer of PDB
                    next if ($file1 =~ /^MASTER|^END/);
                    print $AG_AB_FILE $file1;
                }
    }
        

    
	while (!eof ( $AB_FILE ) )
	{
	    my $file3 = <$AB_FILE>;
	    print $AG_AB_FILE $file3;
            movePDBs ($dir, $destABAG, $pdb_id);
            
	}
                
        movePDBs ($dir, $destABAG, $pdb_id);
        $count++;
        %mapedChains = ();
        
    }
    
    if ( (scalar @antigen) == 2 ) {
        $biAntigen = 1;
    }
    
    return $biAntigen;
       
}



# ************* make_antibody_complex *****************                
# Description: Assembles the antibody chains i.e light and heavy and writes
#              them into a file making antibody complex
# Inputs: A reference to an array with antibody chain labels, PDB code
# Outputs: 
# Subroutine call/Testing: 
# Date: 26 June 2014                                                         
# Author: Saba

sub make_antibody_complex
{
#    my $count = 1;
    my ( $antibody_key_ref, $pdb_name, $hapten, $count, $numbering, $pdb_path,
         $chainsHashRef) = @_;
    my %mapedChains;
    my %chainsHash = %{$chainsHashRef};
    my @LightChains = @{$chainsHash{"Light"}};
    my @HeavyChains = @{$chainsHash{"Heavy"}};
    my ($L, $H);
        
    foreach my $antibody( @$antibody_key_ref )
    {
	my $antibody_LH; 

        if ($hapten)
	{
	    $antibody_LH = $antibody.".pdb";
	}
	else
	{
	    $antibody_LH = $antibody."_num.pdb";
	}
        
        my $outputFile = "$$pdb_name"."_"."$count".".pdb";
        open (my $ABOUT_FILE, '>>', $outputFile) or die "Can not open file\n";
        open (my $ABIN_FILE, '<', $antibody_LH) or die "Can not open file\n";

        my ($L, $H) = split ("", $antibody);
        if ( grep (/$L/, @LightChains) ){
            $mapedChains{L} = $L;
        }
        elsif ( grep (/$L/, @HeavyChains)) {
            $mapedChains{H} = $L;
        }
        if ( grep (/$H/, @LightChains)) {
            $mapedChains{L} = $H;
        }
        elsif ( grep (/$H/, @HeavyChains) ) {
            $mapedChains{H} = $H;
        }
        $mapedChains{A} = [];
        
        printHeader($ABOUT_FILE, $numbering, $pdb_path, %mapedChains);
	while (!eof ( $ABIN_FILE ) )
	{
	    my $file2 = <$ABIN_FILE>;
	    print $ABOUT_FILE $file2;
        }
        close $ABOUT_FILE;
        close $ABIN_FILE;
                
	#rename ( $antibody_LH, "$$pdb_name"."_"."$count".".pdb" );
        $count++;
    }

    return $count;
}

# ************* hasHapten *****************                
# Description: This reads the PDB file and looks for Haptens in the antibody 
#              and returns true if hapten is found while false othrwise
# Inputs:      PDB file
# Outputs:     Boolean 
# Subroutine call/Testing: my $hapten = hasHapten ($file_path); 
# Date: 26 April 2015                                                         
# Author: Saba

 sub hasHapten 
{
    my ($file_path) = @_; 
    my $chaintype = $SFPerlVars::chaintype;
    my $hashapten = $SFPerlVars::hashapten; 
    my $hapten;
    my ($ch, @output);
    my %chainInfo;
    
    # To grab the chain IDs of first 2 proteins in PDB file 
#    my $chains;# = `$chaintype $file_path | head -2`;

    my @chains = `idabchain $file_path`;
    foreach my $chain( @chains) {
        my @chainComp = split (" ", $chain);
        $chainInfo{$chainComp[3]} = substr ($chainComp[2], 0, 1);
                   
    }
    my ($L, $H);
    
    if (( exists $chainInfo{Light} and exists $chainInfo{Heavy}) ) {
        $L = $chainInfo{Light};
        $H = $chainInfo{Heavy};
        @output = `$hashapten -l $L, -h $H $file_path`
    }

    if ( exists $chainInfo{Light}) {
        $L = $chainInfo{Light};
        @output = `$hashapten -l $L $file_path`;
    }

    if ( exists $chainInfo{Heavy}) {
        $H = $chainInfo{Heavy};
        @output = `$hashapten -h $H $file_path`;
    }
          
    if (grep (/^HAPTEN/, @output))
	{
	    $hapten = 1;  
	}
    else
	{
	    $hapten = 0; 
	}
    
    return $hapten; 

}


# ************* processHapten *****************                
# Description: Adds HETATMs (Haptens) to numbered antibodies which are 
#              classified as Hapten containing antibodies
# Inputs:      PDB file, an array refernce with chain IDs of antibody pair
# Outputs:
# Subroutine call/Testing: processHapten($file_path, $hash_keysRef);
# Date: 26 April 2015                                                         
# Author: Saba
sub processHapten
{
    my ($file_path, $antibodyPairsRef) =@_;
    my @antibodyPairs = @ {$antibodyPairsRef};
    
    foreach my $antibodyPair (@antibodyPairs)
    {
	my $antibodyFile = $antibodyPair."_num.pdb";
	my $outputFile = $antibodyPair.".pdb"; 
	`pdbaddhet $file_path $antibodyFile $outputFile`; 
    }
}

sub getchainIDs
    {
        my ( $pdb_path ) = @_;
       my @chains = `idabchain $pdb_path | awk '{print \$3}';`;
       # my @chainIDs = map {substr ($_, 0, 1)} @chains;
       # return @chainIDs;               
    }


sub movePDBs
        {
            my ($dir, $dest, $pdb_id) = @_;
            opendir ( my $DIR, $dir ) or
                die "Can not open directory $dir";
            foreach  my $fl ( grep m/$pdb_id/, readdir ( $DIR ) )
                          {
                              my $from = $dir."/".$fl;
                              move $from, $dest;
                          }
        }

sub getChainTypeWithChainIDs
  {
      my (@chain_info) = @_;
                
    my %myhash = ();
    $myhash{'Antigen'} = [];
    $myhash{'Light'} = [];
    $myhash{'Heavy'} = [];
    my ( @heavy, @light, @antigen, %heavy_light_pair, %ag_hash );

    foreach my $line ( @chain_info )
    {
	#Chain 1, A: Antigen
	my @a=split(' ', $line);
    	my $chainid = $a[2];
	($chainid) = $chainid =~ m/([A-Za-z0-9]{1})/;
      	my $label = $a[3];

	if ( $label eq "Antigen" )
	{
	    push ( @{ $myhash{'Antigen'}}, $chainid );
	    @antigen = @{ $myhash{'Antigen'} };
 	}
	elsif ( $label eq "Light" )
	{
	    push ( @{ $myhash{'Light'}}, $chainid );
	    @light = @{ $myhash{'Light'} };
	}
	elsif ( $label eq "Heavy" )
	{
	    push ( @{ $myhash{'Heavy'}}, $chainid );
	    @heavy = @{ $myhash{'Heavy'} };
	}
    
    }
      return (%myhash);    
  }
            
sub getSingleChainAntibodyAntigenComplex{
    my ( $pdb_id, $destLgAG, $destLg, $numbering, $pdb_path,
         $chainsHashRef, %complex_hash ) = @_;
    my $dir = '.';
    my $count = 1;
    my %chainsHash = %{$chainsHashRef};
    
    my @LightOnlychains = grep {$complex_hash{$_} eq "NULL" } keys %complex_hash;

       
    $count = getsingleChainAntibody(\@LightOnlychains, \$pdb_id, "", $count,
                                   $numbering, $pdb_path, $chainsHashRef);
    movePDBs ($dir, $destLg, $pdb_id);
    my %mapedChains;
    
    foreach my $abChain ( keys %complex_hash )
    {
	my $antibodyChain = $abChain."_num.pdb"; 
	my ($antigenRef, @antigen);
        
        next if ( $complex_hash{ $abChain } eq  "NULL" );
                   
        if ( $complex_hash{ $abChain } ne  "NULL" )
	{
	    $antigenRef = $complex_hash{$abChain};
            @antigen = @{$antigenRef};
        } 

	open ( my $AB_FILE, '<', "$dir/$antibodyChain" ) or 
	    die "Could not open file $antibodyChain";
	open ( my $AG_AB_FILE, '>>',"$pdb_id"."_"."$count".".pdb" ) or
	    die "Can not write complex";

        my @LightChains = @{$chainsHash{"Light"}};
        my @HeavyChains = @{$chainsHash{"Heavy"}};
        my @AntigenChains = @{$chainsHash{"Antigen"}};

        my $Ag = $complex_hash{$abChain};
        
        if ( grep (/$abChain/, @LightChains) ){
            $mapedChains{L} = $abChain;
        }
        elsif ( grep (/$abChain/, @HeavyChains)) {
            $mapedChains{H} = $abChain;
        }
#        if ( grep (/$Ag/, @AntigenChains)) {
 #           $mapedChains{A} = $Ag;
  #      }
        
        my @Ag;
        foreach my $agChain( @antigen) {
            if ( grep (/$agChain/, @antigen)) {
                push (@Ag, $agChain);
            }
        }
        $mapedChains{A} = [@Ag];
        printHeader($AG_AB_FILE, $numbering, $pdb_path, %mapedChains);
        
        foreach  my $agChain(@antigen) {
            
            $agChain = $agChain.".pdb";
            open ( my $AG_FILE, '<', "$dir/$agChain" ) or
                die "Could no open file $agChain";
            
            while ( !eof ( $AG_FILE ) ) 
                {
                    my $file1 = <$AG_FILE>;
                    # To ignore the Footer of PDB
                    next if ($file1 =~ /^MASTER|^END/);
                    print $AG_AB_FILE $file1;
                }
            
        }
        

	while (!eof ( $AB_FILE ) )
	{
	    my $file2 = <$AB_FILE>;
	    print $AG_AB_FILE $file2;
        }
        
        movePDBs ($dir, $destLgAG, $pdb_id);
        $count++;

        
    }

    return $count;

}

sub getResolInfo
    {
        my ($pdbFile) = @_;
        my %resInfo;
        my ($tags, $val);
        
        my @resolInfo = `pdbheader -r -n $pdbFile`;
        foreach my $line (@resolInfo) {
            chomp $line;
            ($tags, $val) = split (/:\s+/, $line);
            $resInfo{$tags} = $val;
            
        }
        return %resInfo;      
        
    }

    
sub printHeader
{
    my ($AG_AB_FILE, $numbering, $pdb_path, %mapedChains ) = @_;
    my %resInfo = getResolInfo($pdb_path);
    my $L = $mapedChains{L};
    my $H = $mapedChains{H};
    my $AgRef = $mapedChains{A};
    
    my @Ag = @{$AgRef};
    my $size = scalar @Ag;
    my ($Ag1, $Ag2);
    
    if ( $size == 1) {
        $Ag1 = $Ag[0];
        $Ag2 = 0;
    }
    elsif ($size > 1) {
        $Ag1 = $Ag[0];
        $Ag2 = $Ag[1];
    }
    
    if ( ( !$Ag1) and (!$Ag2) and ($L) and ($H) )
        {
            print $AG_AB_FILE "REMARK 950 NUMBERING  $numbering\n";
            print $AG_AB_FILE "REMARK 950 METHOD     $resInfo{Type}\n";
            print $AG_AB_FILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
            print $AG_AB_FILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
            print $AG_AB_FILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
            print $AG_AB_FILE "REMARK 950 CHAIN L    $L\n";
            print $AG_AB_FILE "REMARK 950 CHAIN H    $H\n";
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -s $pdb_path`;
        }
    elsif ( ($L) and ($H) and ($Ag1) and (!$Ag2))
        {
            print $AG_AB_FILE "REMARK 950 NUMBERING  $numbering\n";
            print $AG_AB_FILE "REMARK 950 METHOD     $resInfo{Type}\n";
            print $AG_AB_FILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
            print $AG_AB_FILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
            print $AG_AB_FILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
            print $AG_AB_FILE "REMARK 950 CHAIN L    $L\n";
            print $AG_AB_FILE "REMARK 950 CHAIN H    $H\n";
            print $AG_AB_FILE "REMARK 950 CHAIN A    $Ag1\n";
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -s $pdb_path`;
        }
    elsif ( ($L) and ($H) and ($Ag1) and ($Ag2))
        {
            print $AG_AB_FILE "REMARK 950 NUMBERING  $numbering\n";
            print $AG_AB_FILE "REMARK 950 METHOD     $resInfo{Type}\n";
            print $AG_AB_FILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
            print $AG_AB_FILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
            print $AG_AB_FILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
            print $AG_AB_FILE "REMARK 950 CHAIN L    $L\n";
            print $AG_AB_FILE "REMARK 950 CHAIN H    $H\n";
            print $AG_AB_FILE "REMARK 950 CHAIN A    $Ag1\n";
            print $AG_AB_FILE "REMARK 950 CHAIN A    $Ag2\n";
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag2 -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag2 -s $pdb_path`;
            
    }


    
    elsif ( ($L) and (!$Ag1) )
        {
            print $AG_AB_FILE "REMARK 950 NUMBERING  $numbering\n";
            print $AG_AB_FILE "REMARK 950 METHOD     $resInfo{Type}\n";
            print $AG_AB_FILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
            print $AG_AB_FILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
            print $AG_AB_FILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
            print $AG_AB_FILE "REMARK 950 CHAIN L    $L\n";
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -s $pdb_path`;
        }
    elsif ( ($L) and ($Ag1) )
        {
            print $AG_AB_FILE "REMARK 950 NUMBERING  $numbering\n";
            print $AG_AB_FILE "REMARK 950 METHOD     $resInfo{Type}\n";
            print $AG_AB_FILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
            print $AG_AB_FILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
            print $AG_AB_FILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
            print $AG_AB_FILE "REMARK 950 CHAIN L    $L\n";
            print $AG_AB_FILE "REMARK 950 CHAIN A    $Ag1\n";
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $L -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -s $pdb_path`;
        }
    elsif ( ($H) and (!$Ag1))
        {
        print $AG_AB_FILE "REMARK 950 NUMBERING  $numbering\n";
        print $AG_AB_FILE "REMARK 950 METHOD     $resInfo{Type}\n";
        print $AG_AB_FILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
        print $AG_AB_FILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
        print $AG_AB_FILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
        print $AG_AB_FILE "REMARK 950 CHAIN H    $H\n";
        print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -m $pdb_path`;
        print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -s $pdb_path`;
    }
    elsif ( ($H) and ($Ag1))
        {
            print $AG_AB_FILE "REMARK 950 NUMBERING  $numbering\n";
            print $AG_AB_FILE "REMARK 950 METHOD     $resInfo{Type}\n";
            print $AG_AB_FILE "REMARK 950 RESOLUTION $resInfo{Resolution}\n";
            print $AG_AB_FILE "REMARK 950 R-FACTOR   $resInfo{'R-Factor'}\n";
            print $AG_AB_FILE "REMARK 950 R-FREE     $resInfo{'R-Free'}\n";
            print $AG_AB_FILE "REMARK 950 CHAIN H    $H\n";
            print $AG_AB_FILE "REMARK 950 CHAIN A    $Ag1\n";
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $H -s $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -m $pdb_path`;
            print $AG_AB_FILE "REMARK 950 ", `pdbheader -c $Ag1 -s $pdb_path`;
        }

}

sub getsingleChainAntibody
    {
#    my $count = 1;
    my ( $antibody_key_ref, $pdb_name, $hapten, $count, $numbering, $pdb_path,
         $chainsHashRef) = @_;
    my %mapedChains;
    my %chainsHash = %{$chainsHashRef};
    my @LightChains = @{$chainsHash{"Light"}};
    my @HeavyChains = @{$chainsHash{"Heavy"}};
            
    foreach my $antibody( @$antibody_key_ref )
    {
	my $antibody_LH; 

        if ($hapten)
	{
	    $antibody_LH = $antibody.".pdb";
	}
	else
	{
	    $antibody_LH = $antibody."_num.pdb";
	}
        
        my $outputFile = "$$pdb_name"."_"."$count".".pdb";
        open (my $ABOUT_FILE, '>>', $outputFile) or die "Can not open file\n";
        open (my $ABIN_FILE, '<', $antibody_LH) or die "Can not open file\n";

        if ( grep (/$antibody/, @LightChains) ){
            $mapedChains{L} = $antibody;
        }
        elsif ( grep (/$antibody/, @HeavyChains)) {
            $mapedChains{H} = $antibody;
        }
        $mapedChains{A} =[];        

        printHeader($ABOUT_FILE, $numbering, $pdb_path, %mapedChains);
	while (!eof ( $ABIN_FILE ) )
            {
                my $file2 = <$ABIN_FILE>;
                print $ABOUT_FILE $file2;
            }
        close $ABOUT_FILE;
        close $ABIN_FILE;
        
	#rename ( $antibody_LH, "$$pdb_name"."_"."$count".".pdb" );
        $count++;
    }
    
    return $count;
}


        
1;
