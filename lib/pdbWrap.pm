package pdbWrap; 
#***********************************************************************      
#                                                                             
#   Perl Module:   pdbWrap.pm                                                 
#   File:          pdbWrap.pm                                                 
#                                                                             
#   Version:       V1.2                                                       
#   Date:          22.04.15                                                   
#   Fucntion:      The perl module is a libarary of sub routines that are  
#                  wrapper for sub routines from pdb.pm
#                  
#   Modules:       pdb.pm
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
use strict; 
use Exporter qw (import);
use Data::Dumper; 
use Carp;
use File::Copy;
use List::MoreUtils qw(uniq);
use pdb qw 
    (split_pdb_to_chains                                                       
     get_pdbchains_contacts                                                    
     pair_heavy_light_antigen                                                  
     get_hash_key
     antibody_assembly                                                        
     antibody_number
     assemble_CDR_antigen 
     antigen_CDR_conts 
     get_complex 
     get_antibody_antigen_complex 
     make_antibody_complex
     getchainIDs
     getChainTypeWithChainIDs
     getSingleChainAntibodyAntigenComplex
);

our @EXPORT_OK = qw 
    (dirOperations
     processAntibody 
     processAntibodyAntigen 
     movePDBs
     processSingleChainAntibody
     processSingleChainAntibodyAntigen
);

# ************ dirOperations **************                              
# Description: Prepares the directory for PDB processing
#              
# Inputs: Dir name and PDB name 
# Outputs: Working directory name
# Subroutine call/Testing: $dir = dirOperations ($process_dir, $pdb_id);
# Date: 26 April 2015                                                         
# Author: Saba
sub dirOperations
{
    my ($process_dir, $pdb_id) = @_; 
    mkdir $process_dir;
    chdir $process_dir;
    mkdir $pdb_id;
    chdir $pdb_id;
    my $dir = $process_dir.'/'.$pdb_id;
    return $dir; 
}

# ************ processAntibody **************                                 
# Description: Wrapper subroutine from pdb.pm for antibody processing.
#              The description of each subroutine is given before its call 
#                                                                              
# Inputs: PDB file, 2 file handles, numbering scheme flag, PDB code
# Outputs:
# Subroutine call/Testing: my ($antigen, $antigen_ID, $hash_keysRef) =
#                  processAntibody($file_path, $LOG, $SUMMARY, $nsch, $pdb_id);
# Date: 26 April 2015                                                         
# Author: Saba
sub processAntibody
{
    my ($file_path, $LOG, $SUMMARY, $nsch, $pdb_id) = @_;
    my $count_kabatnum_error; 
    
    # Split pdb file into chains 
    split_pdb_to_chains ( $file_path );
    print {$LOG} "$pdb_id is splited into chains\n\n";
    
   # Get hash of inter-chain contacts
    my %contact_hash = get_pdbchains_contacts ( $file_path );
    print {$LOG} "Contacts Summary computed\n\n";
    print {$SUMMARY} "Number of contacs in all chains of $pdb_id\n" .
	Dumper (\%contact_hash);
    
   # Pair heavy and light chains with antigen according to inter-chain contacts
    my ( $heavy_light, $antigen, $antigen_ID, $chainsHash ) =
	&pair_heavy_light_antigen ( $file_path );
    print {$LOG} 
    "Heavy/Light contacts has been saved into seperate hashes".
	"(chain label pair (LH) as keys and number of contacts as values)\n\n";
    print {$SUMMARY} "Number of conatcs between light and Heavy Chains\n".
	Dumper ($heavy_light). "\n". 
	"Number of conatcs between antigen and antigen Chains\n" . 
	Dumper ($antigen) . 
	"All Antigens found in $pdb_id\n" . join ("\n", @$antigen_ID). "\n";

    # Get keys of hash (containing light heavy chains contact information)
    my @hash_keys = get_hash_key ( $heavy_light );
    print {$LOG}"Heavy/Light chains ID are saved to an array for assembly\n\n";

    # Assembling correct antibody (Light and Heavy) and writing into one new 
    # file
    antibody_assembly ( @hash_keys );
    print {$LOG} "Heavy/Light chains are assemebled\n\n";


   return $antigen,$antigen_ID, \@hash_keys, $chainsHash;
 
}

# ************ processAntibodyAntigen **************                       
# Description: Wrapper subroutine from pdb.pm for antibody antigen complex 
#              processing. The description of each subroutine is given before
#              its call    
#                                                                             
# Inputs: Array reference with antigen chain IDs, 2 file handles, antibody
#         chain pairs, PDB code           
# Outputs:                                                                    
# Subroutine call/Testing: processAntibodyAntigen 
#           ($antigen, $antigen_ID, $hash_keysRef, $LOG, $SUMMARY, $pdb_id)

# Date: 26 April 2015                                                         
# Author: Saba
sub processAntibodyAntigen
{
    my ($antigen, $antigen_ID, $hash_keysRef, $LOG, $SUMMARY, $pdb_id,
    $destABAG, $destAB, $dir, $chainsHashRef, $file_path, $numbering) = @_; 
    # Get keys of hash (containg antigen inter-chain contact information) 
    my @antigen_keys = get_hash_key ( $antigen );
        
    # Assembles file containing antibody CRD regions with antigen pdb file
    # Writes them into a separate file (CDR defination + antigen)
    # Computes contacts of CDR regions with each antigen
    my ($nonAgCount, %CDR_hash_ref) = assemble_CDR_antigen($hash_keysRef, $antigen_ID);
    print {$SUMMARY} "Anonymous Hash containg contact information of each".
	"antibody with every antigen" .
	Dumper (\%CDR_hash_ref);
    
    # Compute antibody and antigen complex information and putting that into
    # a hash
    my %complex_hash = get_complex(%CDR_hash_ref);
    print {$SUMMARY} "Hash containg Final complex information\n" .
	Dumper (\%complex_hash);
    
    # Forming final complex, antibody (Light and Heavy) with corresponding 
    # antigen and writing individual complexes in a separate file. Like,
    # 1AFV having 2 complexes will be named as 1AFV_1 and 1AFV_2
    
    get_antibody_antigen_complex ( $pdb_id, $destABAG, $destAB, $chainsHashRef, $numbering,
                                   $file_path, %complex_hash );
    print {$LOG} "$pdb_id file complexes has been written seperatley\n";
    return $nonAgCount;
    
}

# ************ movePDBs **************                            
# Description: Moves the final processed complexes into a new dirctory 
# Inputs: current directory, destination directory, PDB code 
# Outputs:
# Date: 26 April 2015                                                 
# Author: Saba
sub movePDBs
{
    my ($dir, $dest, $pdb_id) = @_;
    opendir ( my $DIR, $dir ) or 
	die "Can not open directory $dir";
    foreach my $fl ( grep m/$pdb_id/, readdir ( $DIR ) )
    {
	my $from = $dir."/".$fl;
	move $from, $dest;
    }
}

sub processSingleChainAntibody
    {
        my ($file_path, $LOG, $SUMMARY, $pdb_id) = @_;
        split_pdb_to_chains ( $file_path );
        print {$LOG} "$pdb_id is splited into chains\n\n";

        my @chainsInfo = `idabchain $file_path`;
        my %chainType_chainID = getChainTypeWithChainIDs(@chainsInfo);
        
        my @antigenChains = @{$chainType_chainID{Antigen}};
        my @lightChains =  @{$chainType_chainID{Light}};
        my @heavyChains =  @{$chainType_chainID{Heavy}};
        
        if ( (@lightChains) and (!@antigenChains) )  {
            return (\%chainType_chainID, \@lightChains);
        }
         elsif ( (@lightChains) and (@antigenChains)) {
            return (\%chainType_chainID, \@lightChains, \@antigenChains);
        }
        elsif ( (@heavyChains) and (!@antigenChains) ){
            return (\%chainType_chainID, \@heavyChains);
        }
        elsif ( (@heavyChains) and (@antigenChains)) {
            return (\%chainType_chainID, \@heavyChains, \@antigenChains);
        }
                
  #       my @Abchains = getchainIDs($file_path);

   #     return (\@Abchains);
    }
    
sub processSingleChainAntibodyAntigen
        {
            my ($ABChainsRef, $AGChainsRef, $pdb_id, $destLgAG, $destLg, $LOG,
                $SUMMARY, $numbering, $pdb_path, $chainsHashRef) = @_;
            my %abchainAg = assemble_CDR_antigen($ABChainsRef, $AGChainsRef);
            print {$SUMMARY} "Anonymous Hash containg contact information of each".
                "antibody with every antigen" .
                    Dumper (\%abchainAg);
            my %lightAntigenPair = get_complex (%abchainAg );
            
            print {$SUMMARY} "Hash containg Final complex information\n" .
                Dumper (\%lightAntigenPair);
            
         
            getSingleChainAntibodyAntigenComplex($pdb_id, $destLgAG, $destLg, $numbering,
                                                 $pdb_path, $chainsHashRef, %lightAntigenPair);
            print {$LOG} "$pdb_id file complexes has been written seperatley\n";
        }
        

1; 
