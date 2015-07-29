#!/acrm/usr/local/bin/perl -s
use strict;
use IO::CaptureOutput qw(capture qxx qxy);
use Carp;
#use warnings;
use Data::Dumper;
use File::Copy;
use Cwd;
use lib "./scripts";
use pdb qw(
    get_pdb_path                                                              
    get_file_data                                                             
    get_pdbcode_list_xml                                                      
    check_chain                                                               
    split_pdb_to_chains                                                       
    get_pdbchains_contacts                                                    
    pair_heavy_light_antigen                                                  
    get_hash_key
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
);
use pdbWrap qw (dirOperations processAntibody processAntibodyAntigen movePDBs);

my $Usage = <<'EOF';
Usage:

    program_name [-k -c -a] [input_file]
 
Options:

numbering scheme

    k : number the antibodies by Kabat numbering scheme
    c : number the antibodies by Chothia numbering scheme
    a : number the antibodies by Martin numbering scheme
    
EOF

if ( defined ( $::help ) ) 
{
    print "$Usage\n";
    exit 0;
}
my ($Antibodies, $Final_Complex);

# Initial numbering scheme variable depending upon user choice on command line
my $nsch;

if ( $::k )
{
    $nsch = '-k';
    $Antibodies  = "Antibody_Kabat";
    $Final_Complex = "Complex_Kabat";
}
elsif ( $::c )
{
    $nsch = '-c';
    $Antibodies  = "Antibody_Chothia";
    $Final_Complex = "Complex_Chothia";
}

elsif ( $::a )
{
    $nsch = '-a';
    $Antibodies  = "Antibody_Martin";
    $Final_Complex = "Complex_Martin";
}

# Define and create process directory 
my $master_dir = getcwd;
my $process_dir = $master_dir.'/'.'Dataprep' . $$ . $nsch; # $$ = process id
mkdir $Antibodies;
mkdir $Final_Complex;

my $input_file = "@ARGV";

# Here are 2 functions,they can be used depending upon type of input file 
# i.e, .xml or .txt
# Use function xml_to_txt, if input file is .xml file
# Use get_file_data, if input file is .txt file

#my @all_pdbs = get_pdbcode_list_xml ( $input_file ) or 
#    die "Can not open $input_file\n";
my @all_pdbs = get_file_data ( $input_file ) or 
    die "Can not open $input_file\n"; 

# Variable Declaration
my $dir;
my $count_superseded = 0;
my $count_kabatnum_error = 0;
my $count_pdb_complex = 0;
my $count_antibody = 0;
my $count_CDR_error = 0;
my $count_hapten = 0; 
my @Only_antibody;
my @antigen_antibody_complex;
my @superseded;
my @kabat_failed;
my @CDR_failed;
my @haptenComplex; 
my %contact_hash;

my ( $MASTER_LOG, $LOG, $SUMMARY );

my (@LgAnt, @HvAnt, @Lg, @Hv, @Ant);
my $LgAnt = 0;
my $HvAnt = 0;
my $Lg = 0;
my $Hv = 0;
my $Ant = 0;

my($light_c, $heavy_c, $antigen_c);
# Open Master file for recording information about Superceded, failed
# Kabat numbering and details of processed PDBs
open ( $MASTER_LOG, ">masterlog".$nsch.".log" );


# Reading the list of antibodies
foreach my $pdb_id ( @all_pdbs )
{
    print "Working on $pdb_id\n"; 
    chomp $pdb_id;

    $dir = dirOperations ($process_dir, $pdb_id);
    # Open log files                                                                                                                                                                                                                       
    open ( $LOG, ">$dir/antibody.log" ) or
        die "Error in creating log file";
    open ( $SUMMARY, ">$dir/cont_summary.txt" ) or
        die "Error in creating log file";
    
    my $file_path = get_pdb_path ( $pdb_id );

# Check if file exists, if it doesn't exist (i.e Superceded) then writes the
# protein ID in Master file.  
    if ( !( -e "$file_path" ) ) 
    {
	push ( @superseded, $pdb_id );
	print {$MASTER_LOG} "$pdb_id is found Superseded\n\n";
	$count_superseded++;
	next;
    }

# Check antibody for types of chains either light, heavy or antigen
# This script only processes 2 types of complexes, 
# 1) Antibody/antigen complexes (light, heavy and antigen)
# 2) Antibody complexes (light and heavy)  
# It ignores all other types of complexes e.g. Only antigen, Only Light,
# Only Heavy and Light or Heavy with antigen. However the details about these 
# are present in master log file

( $light_c, $heavy_c, $antigen_c ) = check_chain ( $file_path );
my $flag = 0; 
my $hapten = hasHapten ($file_path);
if ($hapten)
{
    if ( ( $antigen_c eq "" ) and ( $heavy_c eq 'Heavy' ) and
	     ( $light_c eq 'Light' ) )
    {
	$antigen_c = 'Antigen';
	$flag = 1; 
    }
}

if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq "" ) and 
     ( $light_c eq 'Light') )
{
    push ( @LgAnt, $pdb_id );
    $LgAnt++;
    next;
}

if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq 'Heavy' ) and 
     ( $light_c eq "" ) )
{
    push ( @HvAnt, $pdb_id );
    $HvAnt++;
    next;
}

if ( ( $antigen_c eq "" ) and ( $heavy_c eq 'Heavy' ) and 
     ( $light_c eq "" ) )
{
    push ( @Hv, $pdb_id );
    $Hv++;
    next;
}

if ( ( $antigen_c eq "" ) and ( $heavy_c eq "" ) and 
     ( $light_c eq 'Light' ) ) 
{
    push ( @Lg, $pdb_id );
    $Lg++;
    next;
}

if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq "" ) and 
     ( $light_c eq "" ) ) 
{
    push ( @Ant, $pdb_id );
    $Ant++;
    next;
}


# 1st If-check start
# This if block processes the antibodies with only Light and Heavy chains
if ( ( $antigen_c eq "" ) and ( $heavy_c eq 'Heavy' ) and 
     ( $light_c eq 'Light' ) )
{
    my ($antigen, $antigen_ID, $hash_keysRef) = 
	processAntibody($file_path, $LOG, $SUMMARY, $nsch, $pdb_id);
    
    # Kabatnum is an external program to number antibodies according to 
    # standard numbering scheme - For some antibodies, it does not work and 
    # gives an error. On Failure, the subroutine, antibody_number return an 
    # error by Croak which is being printed on Master log and protein log file.    
    eval { antibody_number ( $hash_keysRef, $nsch ); 1; };
    if ( $@ )
    {
	foreach my $mLOG ($LOG, $MASTER_LOG) 
	{
	    print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
		$@ . "Program exited\n\n";
	}

	push ( @kabat_failed, $pdb_id );
        $count_kabatnum_error++;
	chdir '..';
	next;
    }
    else 
    {
	print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
	    " numbering scheme\n\n";
    }


    # Forming antibody complexes and renaming them depending on number of 
    # complexes. Like, 1AFV having 2 complexes will be named as 
    # 1AFV_1 and 1AFV_2
    make_antibody_complex ( $hash_keysRef, \$pdb_id, $hapten);
        
    # Moving complexes from individual directory to a separate folder to 
    # store all antibody complexes
    my $dest = $master_dir."/". $Antibodies;
    movePDBs ($dir, $dest, $pdb_id); 
    print {$LOG} "Antibody has been moved to the folder of antibody" . 
	"complexes\n\n";
    
    push ( @Only_antibody, $pdb_id );
    $count_antibody++;

} # 1st If-check ends here


# 2nd If-check start                                                          
# This if block processes the antibodies with Light, Heavy and antigen chains 
if ( ( $antigen_c eq 'Antigen' ) and ( $heavy_c eq 'Heavy' ) and 
     ( $light_c eq 'Light' ) )
{ 
    my $dest = "$master_dir"."/".$Final_Complex;
    my ($antigen, $antigen_ID, $hash_keysRef) =
	processAntibody($file_path, $LOG, $SUMMARY, $nsch, $pdb_id);
    
    # Kabatnum is an external program to number antibodies according to 
    # standard numbering scheme - For some antibodies, it does not work and 
    # gives an error. On Failure, the subroutine, antibody_number return an 
    # error by Croak which is being printed on Master log and protein log file.    
    eval { antibody_number ( $hash_keysRef, $nsch ); 1; };
    if ( $@ )
    {
	foreach my $mLOG ($LOG, $MASTER_LOG) 
	{
	    print {$mLOG} "Kabat Numbering Failed for $pdb_id\n" .
		$@ . "Program exited\n\n";
	}
	
	push ( @kabat_failed, $pdb_id );
        $count_kabatnum_error++;
	chdir '..';
	next;
    }
    else 
    {
	print {$LOG} "Antibody has been re-numbered by Kabat/Chotia" .
	    " numbering scheme\n\n";
    }
    if ( ($hapten) and ($flag))
    {
	processHapten($file_path, $hash_keysRef);
	make_antibody_complex ( $hash_keysRef, \$pdb_id, $hapten);
	movePDBs ($dir, $dest, $pdb_id);
	push (@haptenComplex, $pdb_id); 
	$count_hapten++;
	next; 
    }

    # Extract CDRs regions according to Kabat CDR definition, eval function 
    # check if CDR extraction get failed Writes CDR definitions into a 
    # seperate file
    eval { extract_CDR ( $hash_keysRef ); 1; };
    if ( $@ )
    {
	foreach my $mLOG ( $LOG, $MASTER_LOG )
	{
	    print {$mLOG} "Can'nt extract CDRs for: $pdb_id\n" . $@ .
		"Program exited\n\n";
	}
	push ( @CDR_failed, $pdb_id );
        $count_CDR_error++;
	chdir '..';
	next;
    }
    else
    {
	print {$LOG} "CDRs for $pdb_id has been extracted from each light/Heavy".
	    "chain\n\n";
    }
    
    processAntibodyAntigen ($antigen, $antigen_ID, $hash_keysRef, $LOG, $SUMMARY, $pdb_id);

    # Moving complexes from individual directory to a separate folder to 
    # store all antibody/antigen complexes 
    #my $dest = "$master_dir"."/".$Final_Complex;
    movePDBs ($dir, $dest, $pdb_id);
    
    push ( @antigen_antibody_complex, $pdb_id );
    $count_pdb_complex++;
} # 2nd If-check ends here

chdir '..';
#last; 
} # Main For loop ends here

print "$count_pdb_complex proteins has been successfully processed " .
    "for antibody/antigen complexes\n";
print "$count_antibody proteins were found as only antibody complexes\n";
print "$LgAnt proteins were found as Light-Antigen complexes\n";
print "$HvAnt proteins were found as Heavy-Antigen complexes\n";
print "$Lg proteins were found as Light chains\n";
print "$Hv proteins were found as Heavy chains\n";
print "$Ant proteins were found as Antigen chains (Fc Fragments)\n";
print "$count_kabatnum_error proteins failed to number by program Kabatnum\n";
print "$count_superseded proteins were found as Superseded\n";
print "$count_CDR_error proteins were failed for CDR extraction\n";
print "$count_hapten antibodies are complexed with Hapten\n"; 
#######################

chdir '..';

# open files  for recording files
open (my $ABAGCOMPLEX, '>Antigen_Antibody_complexes.list') or
die "Can not open file/n";

open (my $FREE_ANTIBODY, '>Antibody_complexes.list') or
die "Can not open file/n";

open (my $LG_AG, '>Light_Antigen_complexes.list') or
die "Can not open file/n";

open (my $HV_AG, '>Heavy_Antigen_complexes.list') or
die "Can not open file/n";

open (my $LG, '>Light_Only.list') or
die "Can not open file/n";

open (my $HV, '>Heavy_Only.list') or
die "Can not open file/n";

open (my $FC, '>FC_Fragments.list') or
die "Can not open file/n";

open (my $FAILED, '>Kabat_Failed.list') or
die "Can not open file/n";

open (my $SUPERSEDED, '>Superseded.list') or
die "Can not open file/n";

open (my $CDR_ERROR, '>CDR_Error.list') or
die "Can not open file/n";

open (my $HAPTEN, '>Hapten_Complex.list') or 
    die "Can not open file\n"; 

print {$MASTER_LOG} ">$count_pdb_complex proteins has been successfully " . 
    "processed for antibody/antigen complexes\n";
print {$MASTER_LOG} "@antigen_antibody_complex\n";
print {$ABAGCOMPLEX} join("\n", @antigen_antibody_complex);
print {$MASTER_LOG} ">$count_antibody proteins were found as only antibody " .
    "complaxes\n";
print {$MASTER_LOG} "@Only_antibody\n";
print {$FREE_ANTIBODY} join ("\n", @Only_antibody), "\n";

print {$MASTER_LOG} ">$LgAnt proteins were found as Light-Antigen complexes\n";
print {$MASTER_LOG} "@LgAnt\n";
print {$LG_AG} join ("\n", @LgAnt), "\n";

print {$MASTER_LOG} ">$HvAnt proteins were found as Heavy-Antigen complexes\n";
print {$MASTER_LOG} "@HvAnt\n";
print {$HV_AG} join("\n", @HvAnt), "\n";

print {$MASTER_LOG} ">$Lg proteins were found as Light chains\n";
print {$MASTER_LOG} "@Lg\n";
print {$LG} join ("\n", @Lg), "\n";

print {$MASTER_LOG} ">$Hv proteins were found as Heavy chains\n";
print {$MASTER_LOG} "@Hv\n";
print {$HV} join ("\n", @Hv), "\n";

print {$MASTER_LOG} ">$Ant proteins were found as Antigen chains (Fc Fragments)\n";
print {$MASTER_LOG} "@Ant\n";
print {$FC} join ("\n", @Ant), "\n";

print {$MASTER_LOG} ">$count_kabatnum_error proteins failed to number by" .
    "program Kabatnum\n";
print {$MASTER_LOG} "@kabat_failed\n"; 
print {$FAILED} join ("\n", @kabat_failed), "\n"; 

print {$MASTER_LOG} ">$count_superseded proteins were found as superseded\n";
print {$MASTER_LOG} "@superseded\n";
print {$SUPERSEDED} join ("\n", @superseded), "\n";

print {$MASTER_LOG} ">$count_CDR_error proteins were failed for CDR extraction".
    "\n";
print {$MASTER_LOG} "@CDR_failed\n";
print {$CDR_ERROR} join ("\n", @CDR_failed), "\n";

print {$MASTER_LOG} "@haptenComplex";
print {$HAPTEN} join ("\n", @haptenComplex), "\n";

my $total = $count_pdb_complex+$count_antibody, "\n";

print {$MASTER_LOG} "Complexes=$count_pdb_complex\n";
print {$MASTER_LOG} "Free_Antibody=$count_antibody\n";
print {$MASTER_LOG} "Complete_Dataset=$total\n";

print {$MASTER_LOG} "Light-Antigen=$LgAnt\n";
print {$MASTER_LOG} "Heavy-Antigen=$HvAnt\n";
print {$MASTER_LOG} "Bens-Jones=$Lg\n";
print {$MASTER_LOG} "Camelids=$Hv\n";
print {$MASTER_LOG} "Fc=$Ant\n";
print {$MASTER_LOG} "Failed=$count_kabatnum_error\n";
print {$MASTER_LOG} "Superseded=$count_superseded\n";
print {$MASTER_LOG} "Hapten=$count_hapten\n";
print "See masterlog.log for details...\n";

=head1 NAME

get_antibody_complex.pl 

=head1 SYNOPSIS

To use this program, type on command line, 
./complex_antibody.pl [-k -c -a] input file

=head1 DESCRIPTION

 Perl Module for generating separate antibody complexes, present in a PDB file
 This prgram has an option of 
 numbering antibodies by 3 different numbering schemes namely Kabat, Chothia 
 and Martin.

=head1 AUTHOR

Saba Ferdous (ucbterd@acrm19)

=head1 COPYRIGHT AND LICENSE
Copyright (C) 2014 by Saba Ferdous

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS
None reported... yet.

=cut
