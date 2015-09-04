# This perl module has the running subroutines for the idea of obtaining nearest 
# point to the MEDIAN of the line and then forward trace and back trace.  

package epitope;

use strict;
#use warnings;
 
use Exporter qw(import);
my $pymol = "/usr/bin/pymol"; 
our @EXPORT_OK = qw(getRegionRange getFragmentResidue writepymolScript 
getDistanceToStraightLine getChainLabelAndType getCACoordsForRegion getXYZ
calculateDistanceWithLine getBestFitLineCoords makePeptideFigure
SecondaryStrAssignmentToPeptides printArrays getLeastDistanceCAlpha 
getEuclideanDistance SumOfDistance getSide getPointfromDistance 
getPercentAlpha makePeptideFigure getVectors getVectorAngle
getVecPointByArrayofPoints getMidPoints absoluteDistancePoints
getLeastDistanceCAlpha getNearestCAToMidPoint getSecondclosestpoint 
getSide getPointfromDistance getReferencepoints getAllCADistanceWithLine
getPointsForCAonLine getPointToPointDistanceOnLine averageDeviation 
reverseRegLineDir checkDeviationOnEnds getResDist getResSeq getContacts);
use Math::Vec qw(NewVec);
use SFPerlVars;
use List::Util qw(reduce);
use File::Basename;
use Data::Dumper;
use Math::Vector::Real;
use Math::Trig;

# Description: A subroutin to find the start and end element of an array
# Inputs: A string of epitope residues (regions) seperated by commas
# Outputs: start and end element of each of epitope (array of anonymous arrays)
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub getRegionRange
{
    my ($regions) = @_;
    my (@splittedRegion, @rangeRegion, @regions);
    # Each region is separated by comma, For example: 
    # 328 329 330 331 332,341 342 343 344,366 368 369 370 372,400 401 402 403
   
    @regions = split (/,/, $regions);

    foreach my $region(@regions)
    {
        @splittedRegion = split (" ", $region);
        my $lastElement = pop @splittedRegion; # To get last element
	# Region Range = 328-332; start to end
        push ( @rangeRegion, [$splittedRegion[0], $lastElement ] )
    }

    # Array of anonymous arrays i.e containing addresses of arrays        
    return @rangeRegion;            

}

# Description: This subroutine finds the residues of fragments
# Inputs: A string of epitope residues (fragments) seperated by commas
# Outputs: An array containg all residues forming fragments
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub getFragmentResidue
{
    my ($fragments) = @_;
    my (@splittedFragments, @fragments);

    # Each fragment is seperated by comma and further 
    # each residue in fragment is seperated by space 
    # For Example: 90, 102 103

    @fragments = split (/,/, $fragments);
    
    foreach my $fragment (@fragments)                                        
    {   
	chomp $fragment;
	push (@splittedFragments, split (" ", $fragment) );                  
    }

    return  @splittedFragments;
}

# Description: This subroutine generates a figure in pymol where epitope 
#              regions are highlighted in red colour while fragments in green 
#              color
# Inputs: 2 array references; 1) containg range of all the regions present in 
#         epitope. 2) containing all the residues forming fragments.3) Antigen
#         chain label. 4) PDB (aligned) file bringing epitope region on front
#         5) original PDB file (Antibody- Antigen complex)    
# Outputs: A figure in png format
# Subroutine call/Testing: 
# Date: 03 NOv 2014                                 
# Author: Saba

sub writepymolScript
{
    my ($rangeRegionRef, $residueFragmentRef, 
	$antigenChainLabel, $alignedPdb, $pdbFile) = @_;
    my $pdbFileName = basename($pdbFile, ".pdb"); 
    
    # Writing pymol script to provide as input to pymol program
    open (my $PYMOL, '>', "pymolscript.pml") or die
	"Can not open file"; # opening perl script file   
    
    print {$PYMOL} <<__EOF;
load $alignedPdb                     
bg_color white                                                
turn x, 90                                                    
turn y, 90                                                    
Turn x, 90                                                    
turn y, 90                                                    
show cartoon                                                  
hide lines                                                    
select light, chain L                                         
remove light                                                  
select heavy, chain H                                         
remove heavy                                                  
color cyan, chain $antigenChainLabel
__EOF

# Select the given epitope range and color that in red
foreach my $record (@{$rangeRegionRef})
{ 
    # Accessing elements from anonymous array
    print {$PYMOL} "select regions, resi $record->[0]-$record->[1]\n";
    print {$PYMOL} "color red, regions\n";
}

# Select and color the fragment residues
foreach my $record (@{$residueFragmentRef})
{
    print {$PYMOL} "select fragments, resi $record\n";
    print {$PYMOL} "color green, fragments\n";
}
    
# Taking the image
print {$PYMOL} "ray 800,600\n";
print {$PYMOL} "png $pdbFileName.png\n";
print {$PYMOL} "quit\n";

# Sending pymol script to program pymol
`$pymol -c pymolscript.pml`;
unlink "pymolscript.pml";
unlink $alignedPdb;
    
}

############################ Distance #####
# Description: 
# Inputs: 
# Outputs: 
# Subroutine call/Testing: 
# Date: 2014                                 
# Author: Saba

sub getDistanceToStraightLine
{
    my ($rangeRegionRef, $PDBAllCoordsRef, $antigenChainLabel, $pdbFile, 
	$SUMRY, $dir, $REGLEN, $CURVEDLEN, $LINEARLEN, $FOLDEDLEN, $SSINFO) = @_;
    my $coords;
    my @coordsCA = ();
    #my $P; 
    my ($P, $P0, $P1, $angle);
    my $count = 1;
    my @AlldistToLine = ();
    my @AllepitopeSS = ();
    my @AllDeviations = ();
    my @AllpointsToMidpoint = ();
    my @epitopeSS;
    my ($CA, $closestCA);
    
    my $pdbFileName = basename($pdbFile, ".pdb");    
    # Loop reads (perform operations on) each region every time
    print {$REGLEN} "$pdbFile: ";
    
    foreach my $record (@{$rangeRegionRef})
    {
	print {$SUMRY} "region_".$count."\n";

	my ($startRes, $endRes) = ($record->[0], $record->[1]);
	# Secondary structure assignment to each region
	@epitopeSS = SecondaryStrAssignmentToPeptides 
	    ($pdbFile, $antigenChainLabel, $startRes, $endRes, $count);
	
	print {$SUMRY} "Secondry sturcture of region\n";
	print {$SUMRY} "@epitopeSS\n";
	

	# Obtains CA atom records for each region
	@coordsCA = getCACoordsForRegion
	    ($PDBAllCoordsRef, $antigenChainLabel, $startRes, $endRes);
	my $sizePep =  scalar @coordsCA;
	print {$SUMRY} "Size of peptide = $sizePep\n";
        print {$REGLEN} "$sizePep ";
        
	print {$SUMRY} "CA coordinates of the region has been obtained\n";
	print {$SUMRY} Dumper (\@coordsCA), "\n";

	# Output file name with region range
	my $regressionOF = "regression".$startRes."_".
	    $endRes.$pdbFileName.".pdb";

	# Get Best Fit line for the peptide of given range
	getBestFitLine($startRes, $endRes, $antigenChainLabel,
		       $pdbFile, $regressionOF); 
	print {$SUMRY} "Best fit line in the region is obtained\n";
	
	# Making figure of region with line of best fit in pymol 
	makePeptideFigure("$regressionOF", $count);

	# Opening this file before sending to subroutine
	# Sending file handle does not require file to be in the same directory
	# where the script is being run from
	open (my $IN, "$regressionOF") or die "Can not open file";
	my @regressionCoords = getBestFitLineCoords($IN); 

	print {$SUMRY} "Co-ordinates of best fit line are obtained\n";
	print {$SUMRY} Dumper (\@regressionCoords), "\n";
	
	my ($P0, $P1, $VReg) = getRegressionVector(@regressionCoords); 
	my $VPep = getPeptideVector(@coordsCA); 
	print {$SUMRY} "Vector is obtained by the two given points\n";
	print {$SUMRY} "VReg = $VReg\n";
	print {$SUMRY} "VPep = $VPep\n";
	
	# Seek function is used to take the filehandle again 
	# on the start of file
	seek ($IN, 0, 0); 

	# Get resSeq to work out for contacts                                                                      
        my @resSeq = getResSeq($IN);

        $angle = getVectorAngle ($VReg, $VPep);
	
	print "The Angle is $angle\n";
	print {$SUMRY} "The Angle is $angle\n";


        # Seek function is used to take the filehandle again                                                      
        # on the start of file                                                                                     
        seek ($IN, 0, 0);

	# To make sure that the 2 vectors, regression line vactor and 
	# peptide vector are in the same direction. If 2 vectors are in same 
	# direction then angle between them would be less then 90 otherwise
	# vectors are in opposite direction. In such case, we need to reverse
	# the direction of regression line. 
	if ($angle > 101)
	{
	    reverseRegLineDir ($IN, $regressionOF);
	    open (my $IN2, $regressionOF) or die "Cannot open\n"; 

	    my @regressionCoords = getBestFitLineCoords($IN2);

	#    print "DEBUGGG: @regressionCoords\n";
#exit;  	    
	    print "The regression line direction has been reversed\n\n\n";
	    print {$SUMRY} "The regression line direction has been reversed\n\n";
	    my $VRRev;
	    ($P0, $P1, $VRRev) = getRegressionVector(@regressionCoords);
	    $angle = getVectorAngle ($VRRev, $VPep);
	    print "New angle with reversed regression line: $angle\n";
	    print {$SUMRY} "New angle with reversed regression line: $angle\n";
	}

	seek ($IN, 0, 0);
	my $midPoint = getMidPoints ($P0, $P1);

	print {$SUMRY} "The mid point on Best fit line has been mapped\n";
	print {$SUMRY} "Mid point = $midPoint\n";
	
	# To find the nearest CA to the midpoint
	my ($CAIndex, $closestCADistance, @MidpointToCAdistances) =
	    getNearestCAToMidPoint ($midPoint, @coordsCA);
	print {$SUMRY} "Absloute distances between CA and midpoint is obtained\n";
	print {$SUMRY} "Distances of all CAs with midpoint are:\n",
	join ("\n",  @MidpointToCAdistances), "\n";
	print {$SUMRY} "Index of closest CA = $CAIndex\n";
	print {$SUMRY} "Distance of closest CA = $closestCADistance\n";

	# To avoid the first and last CA as the nearest CA
	# Because due to fold/bend at the end, these could be closest
	# therefore, in such cases, the second closest point would be taken 
	# as the closest point to CA
	if ( ($CAIndex == $sizePep) or ($CAIndex == 1) ) 
	{
	    ($CAIndex, $closestCADistance) =
		getSecondclosestpoint ($CAIndex, @MidpointToCAdistances);
	    
	# If last point is the closet point then look for second closet point
	# e.g 4DGI_1
	    print {$SUMRY} "Index of SECOND closest CA = $CAIndex\n";
	    print {$SUMRY}"Distance of SECOND closest CA = $closestCADistance\n";
	}
	# To map the ideal/reference points on the regression line
	my  ($ss, %refPoints) =	
	    getReferencepoints(\@coordsCA, $CAIndex, $closestCADistance,
			       $P0, $P1, $regressionOF,  \@epitopeSS, $SUMRY, $SSINFO);

	print {$SUMRY} "Hash contaning all the reference points has been ". 
	    "obtained\n";
	print {$SUMRY} Dumper (\%refPoints);

	# RF = First reference point on line
	# RL = Last reference point on line
	# Here First and last point of regression line could also be used but 
	# I used first and last reference point due to the fact if line is 
	# shorter than the the peptide
	my @RF = split(',', $refPoints{R1});
	my @RL = split(',', $refPoints{"R".$sizePep});  

	# Deviation on Straight Line
	my ($RefpointsVec, @CAdistToLine) = 
	    getAllCADistanceWithLine (\@RF, \@RL, \@coordsCA, \%refPoints);

	print {$SUMRY} "Vector of Best fit line (based on first and last ".
	    "reference point): $RefpointsVec\n";
	print {$SUMRY} "Distances of CA to lines have been obtained\n";
	print {$SUMRY} "CA to Line distance is :", join ("\n", @CAdistToLine),
	"\n";

	# To map the actual CA points on the regression line
	my %CApointsOnLine = 
	    getPointsForCAonLine ($RefpointsVec, \@RF, \@coordsCA, 
				  \@CAdistToLine, \$regressionOF); 
	
	print {$SUMRY} "Hash contaning all the CA actual points has been ".
	    "obtained\n";
	print {$SUMRY} Dumper (\%CApointsOnLine);
	print {$SSINFO} "$pdbFile:";
        # To map the deviation between ideal point and actual CA point
	# RF-CA distance on line provides deviation
	my @IdealActualDeviations = 
	    getPointToPointDistanceOnLine(\%refPoints, \%CApointsOnLine); 

	print {$SUMRY} "ideal and actual point deviations obtained\n";
	print {$SUMRY} join("\n",@IdealActualDeviations), "\n";

	# Sum of all the distances between CAs and reference points
	my $averageDeviation = 
	    averageDeviation(@IdealActualDeviations);
	
	my $pdbFileName = basename($regressionOF, ".pdb");
	
 	my ($maxcontactsLocal, $maxcontactsDistant, $totalContacts) =
	    getContacts ($antigenChainLabel, $pdbFile, $sizePep,
			 $ss, $SUMRY, @resSeq); 
	
	classifyShape($averageDeviation, $sizePep, $pdbFile, $count, $dir, 
		      $SUMRY, $pdbFileName, $regressionOF, $angle, 
		      $maxcontactsLocal, $maxcontactsDistant,
		      $totalContacts, $CURVEDLEN, $LINEARLEN, $FOLDEDLEN, @IdealActualDeviations); 
	
	print {$SUMRY} "The average deviation is: $averageDeviation\n";

	push (@AllepitopeSS, [@epitopeSS]);
	push (@AllDeviations, [@IdealActualDeviations]);
        $count++;
    }
    print {$REGLEN} "\n";
    
    return (\@AllepitopeSS, \@AllDeviations);
}

sub classifyShape
{
    my ($averageDeviation, $sizePep, $pdbFile, $count, $dir, $SUMRY,
        $pdbFileName, $regressionOF, $angle, $maxcontactsLocal, 
	$maxcontactsDistant, $totalContacts, $CURVEDLEN, $LINEARLEN, $FOLDEDLEN, @IdealActualDeviations)
	= @_; 
    my $pepFigureFile = "$pdbFileName.png";
    my $pepPDBFile = "$regressionOF"; 
    
#    print "$pepPDBFile\n"; 
 #   exit; 

    my ($avgDev, $shape); 
    if ($averageDeviation <= 1.0)
    {
	if ($sizePep <= 4)
	{
	    if ($averageDeviation < 0.6)
	    {
		print "The peptide region_".$count," of $pdbFile is LINEAR\n";
		print {$SUMRY} "The peptide region_".$count," of $pdbFile is".
		    " LINEAR\n";
		`cp $pepFigureFile $pepPDBFile $dir/LINEAR`;
                print {$LINEARLEN} "$pdbFile:$sizePep\n";
                
	    }
	    else
	    {
		if ( ($maxcontactsDistant >= 2) or ($maxcontactsLocal >= 3 ) or
		     ($totalContacts >=3) )
		{
		    $shape = "folded";
		}
		else
		{
		    $shape = "curved";
		}

#		$shape = checkAngle($angle);
		print {$SUMRY} "The peptide region_".$count," of $pdbFile ".
		    "is $shape\n";
		print "The peptide region_".$count," of $pdbFile is $shape\n";
		if ($shape eq "curved")
		{
		    `cp $pepFigureFile $pepPDBFile $dir/CURVED`;
                    print {$CURVEDLEN} "$pdbFile:$sizePep\n";
                    
		}
		elsif ($shape eq "folded")
		{
		    `cp $pepFigureFile $pepPDBFile $dir/FOLDED`;
                    print {$FOLDEDLEN} "$pdbFile:$sizePep\n";
		}
	    }
	}
	
	else
	{	    
	    print "The peptide region_".$count," of $pdbFile is LINEAR\n";
	    print {$SUMRY} "The peptide region_".$count," of $pdbFile is ".
		"LINEAR\n";
	    `cp $pepFigureFile $pepPDBFile $dir/LINEAR`;
            print {$LINEARLEN} "$pdbFile:$sizePep\n";
        }
	
    }
    elsif (($averageDeviation > 1.0) and ($averageDeviation < 2.5))
    {
	if ($sizePep >= 6)
	{
	    ($avgDev) = 
		checkDeviationOnEnds (@IdealActualDeviations);
	    print "The average deviation of hooked peptide is = $avgDev\n"; 
	    print {$SUMRY} "The average deviation of hooked peptide is =".
		" $avgDev\n";
	    
	    if ($avgDev <= 1.0)
	    {
		if ($avgDev >= 0.5)
		{
		    if ( ($IdealActualDeviations[0] > 1.0 ) and
			 ($IdealActualDeviations[$#IdealActualDeviations] 
			  > 1.0) )  
		    {
			print "The peptide region_".$count," of $pdbFile is".
			    " CURVED\n";
			print {$SUMRY} "The peptide region_".$count," of ".
			    "$pdbFile is CURVED\n";
			`cp $pepFigureFile $pepPDBFile $dir/CURVED`;
                        print {$CURVEDLEN} "$pdbFile:$sizePep\n";
                        
                    }
		    else 
		    {
			print "The peptide region_".$count," of $pdbFile is ".
			    "LINEAR\n";
			print {$SUMRY} "The peptide region_".$count," of ".
			    "$pdbFile is LINEAR\n";
			`cp $pepFigureFile $pepPDBFile $dir/LINEAR`;
                        print {$LINEARLEN} "$pdbFile:$sizePep\n";
                    }

		}
		else
		{
		    print "The peptide region_".$count," of $pdbFile is ".
			"LINEAR\n";
		    print {$SUMRY} "The peptide region_".$count," of $pdbFile".
			" is LINEAR\n";
		    `cp $pepFigureFile $pepPDBFile $dir/LINEAR`;
                    print {$LINEARLEN} "$pdbFile:$sizePep\n";
                }
	    }
	    
	    elsif ($avgDev > 1.0) 
	    {

		if ( ($maxcontactsDistant >= 2) or ($maxcontactsLocal >= 3 ) or
		     ($totalContacts >=3) )
		{
		    $shape = "folded";
		}
		else
		{
		    $shape = "curved";
		}

		#$shape = checkAngle($angle);
		print {$SUMRY} "The peptide region_".$count," of $pdbFile is ".
		    "$shape\n";  
		print "The peptide region_".$count," of $pdbFile is $shape\n";
		if ($shape eq "curved")
		{
		    `cp $pepFigureFile $pepPDBFile $dir/CURVED`;
                    print {$CURVEDLEN} "$pdbFile:$sizePep\n";
                }
		elsif ($shape eq "folded")
		{
		    `cp $pepFigureFile $pepPDBFile $dir/FOLDED`;
                    print {$FOLDEDLEN} "$pdbFile:$sizePep\n";
                }
	    }
	}
	else
	{
#	    $shape = checkAngle($angle);
	    if ( ($maxcontactsDistant >= 2) or ($maxcontactsLocal >= 3 ) or
		 ($totalContacts >=3) )
	    {
		$shape = "folded";
	    }
	    else
	    {
		$shape = "curved";
	    }
	    
	    print {$SUMRY} "The peptide region_".$count," of $pdbFile is ".
		"$shape\n";
	    print "The peptide region_".$count," of $pdbFile is $shape\n";
	    if ($shape eq "curved")
	    {
		`cp $pepFigureFile $pepPDBFile $dir/CURVED`;
                print {$CURVEDLEN} "$pdbFile:$sizePep\n";
            }
	    elsif ($shape eq "folded")
	    {
		`cp $pepFigureFile $pepPDBFile $dir/FOLDED`;
                print {$FOLDEDLEN} "$pdbFile:$sizePep\n";
            }
	    
	}
    }
    
    elsif ($averageDeviation > 2.5)
    {
	if ( ($maxcontactsDistant >= 2) or ($maxcontactsLocal >= 3 ) or 
	     ($totalContacts >=3) )
	{
	    $shape = "folded";
	}
	else
	{
	    $shape = "curved";
	}
	
	print {$SUMRY} "The peptide region_".$count," of $pdbFile is $shape\n";
	print "The peptide region_".$count," of $pdbFile is $shape\n";        
	if ($shape eq "curved")                                               
	{   
	    `cp $pepFigureFile $pepPDBFile $dir/CURVED`;                     
            print {$CURVEDLEN} "$pdbFile:$sizePep\n";
        }                                                                     
	elsif ($shape eq "folded")                                            
	{                                                                     
	    `cp $pepFigureFile $pepPDBFile $dir/FOLDED`;   
            print {$FOLDEDLEN} "$pdbFile:$sizePep\n";
        }      
    }
    
}


sub getBestFitLine
{
    my ($startRes, $endRes, $antigenChainLabel,
        $pdbFile, $regressionOF) = @_;
    my $pdbline = $SFPerlVars::pdbline;
    
    my $startRes = $antigenChainLabel.".".$startRes; # e.g: A23               
    my $endRes = $antigenChainLabel.".".$endRes; # e.g: A28                   
    
    # Using program pdbline to obtain line of best fit for each region        
    `$pdbline -r LIN -a O $startRes $endRes ../$pdbFile $regressionOF`;
}

sub getRegressionVector
{
    my (@regressionCoords) = @_;
    my @startReg =  @{$regressionCoords[0]};
    my @endReg = @{$regressionCoords[$#regressionCoords]};
    my $VReg = getVectors (\@startReg, \@endReg);
    
    my $P0 = getVecPointByArrayofPoints (@startReg);
    my $P1 = getVecPointByArrayofPoints (@endReg);
    
    return ($P0, $P1, $VReg);
}

sub getPeptideVector
{
    my (@coordsCA) = @_;
    my @startPep = @{$coordsCA[0]};
    my @endPep = @{$coordsCA[$#coordsCA]};
    my $VPep = getVectors (\@startPep, \@endPep);
    return $VPep;
}


# Description: Finds the chain type and chain label of the first chain in PDB,
#              In my Complexes, first chain is always antigen
# Inputs: PDB file (Antibody-Antigen Complex)
# Outputs: Antigen chain type (Protein/DNA/RNA), Antigen chain label
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub getChainLabelAndType
{
    my ($pdbFile) = @_;
    my $chaintype = $SFPerlVars::chaintype;
    
    # To obtain the antigen chain label and chain type (N Protein)         
    my ($antigenChainLabel, $antigenChainType) =
	split(" ", `$chaintype $pdbFile | head -1`);
    return $antigenChainLabel, $antigenChainType;
}

# Description: Obtains CA co-ordinates for a given range of an epitope 
# Inputs: PDB file coordinates
# Outputs: CA co-ordinates for given range of PDB (epitope)
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub getCACoordsForRegion
{
    my ($PDBAllCoordsRef, $antigenChainLabel, $startRes, $endRes) = @_;
    my (%parsedPDB, @CACoords, $coords, $x, $y, $z);
 
   # map function works like a foreach loop and returns an array
    my @atomHrefs = map { {parseATOMLine($_)} } @{$PDBAllCoordsRef};
    # each line from $PDBAllCoordsRef is returned as a hash reference and 
    # @atomHrefs contains references of each hash 

    my $startResID = $antigenChainLabel . $startRes;
    my $endResID   = $antigenChainLabel . $endRes;

   # print "DEBUGG: $startResID\n";
  #  print "DEBUGG: $endResID\n";
    
    addResIDsToAtoms(@atomHrefs);
    addResidueIndexToAtoms(@atomHrefs);
    
#    print "resid " . $_->{resID} . " index " . $_->{residueIndex} . "\n"
#	foreach @atomHrefs;
   #exit; 
    my $startIndex = findResidueIndexForResID($startResID, @atomHrefs);
    my $endIndex   = findResidueIndexForResID($endResID, @atomHrefs);
    #print "DEBUGG: $startIndex\n";
    #print "DEBUGG: $endIndex\n";
    #exit; 
    my @regionCAHrefs = 
	map {[$_->{x}, $_->{y}, $_->{z}]} 
    grep { isCAInRange($_, $startIndex, $endIndex) } @atomHrefs;
    
    return @regionCAHrefs;  
}


sub addResIDsToAtoms {
    my (@atomHrefs) = @_;
    #my %atomHrefs;    
    foreach my $atomHref (@atomHrefs) {
        $atomHref->{resID} =
            $atomHref->{chainID}.$atomHref->{resSeq}.$atomHref->{iCode};
        # Concatenate chainid + resSeq + iCode                               
    }
}

sub addResidueIndexToAtoms {
    my @atomHrefs = @_;

    my $prevResID = "";
    my $residueIndex = 0;
    
    foreach my $atomHref (@atomHrefs){
        my $resID = $atomHref->{resID};
        if($resID ne $prevResID){
            ++$residueIndex;
            $prevResID = $resID;

        }
	$atomHref->{residueIndex} = $residueIndex;
    }
}

sub findResidueIndexForResID {
    my $resID = shift;
    my @atomHrefs = @_;
    
    foreach my $atom (@atomHrefs){
        if($atom->{resID} eq $resID){
            return $atom->{residueIndex};
        }
    }
    # If reached here, no atom matching resID                                
}

sub isCAInRange {
    my ($atomHref, $rangeStartIndex, $rangeEndIndex) = @_;

    return  $atomHref->{name} eq 'CA' && 
	indexInRange($atomHref->{residueIndex}, $rangeStartIndex,
		    $rangeEndIndex);
#&& $atomHref->{chainID}


}

sub indexInRange {
    my ($index, $startRangeIndex, $endRangeIndex) = @_;

    return $index >= $startRangeIndex && $index <= $endRangeIndex;
}


# Description: Assigns x y and z co-ordinates to local variables 
# Inputs: Co-ordinates of a Vector point 
# Outputs: 3 Variables containing x y and z co-ordinates
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba
sub getXYZ
{
    my ($startP) = @_;
    my ($x, $y, $z) = split (/\s+/, $startP);
    return $x, $y, $z;
}

# Description: Finds distance of all CA of given epitope with straight line
# Inputs: Start and End points of a line and CAs of epitope
# Outputs: 1) 2 vectors; i) Stright line - $VL ii) Point on epitope and first 
#         point of straight line  - $W 2) array of distances for epitope CA 
#         and straight line 
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub calculateDistanceWithLine
{
    my ($P0, $VL, @coordsCA) = @_;
    my ($P, $W, @CAdistToLine, $distance);
    
    # This for loop reads all CAs including first and last 
    for (my $i = 0; $i<= (scalar @coordsCA)-1 ; $i++)
    {
	# To get co-ordinates of first point (CA) of epitope 
	my @xyzCoords = @{$coordsCA[$i]}; 
	$P = getVecPointByArrayofPoints(@xyzCoords); 
	$W = $P - $P0; # Vector betweem CA (P) and first point on line (P0) 
	
	$distance = pointToLineDistance ($VL, $W); 
	push (@CAdistToLine, $distance);
    }
    return (@CAdistToLine);
}

sub pointToLineDistance
{
    my ($VL, $W) = @_;
    my ($product, $magnitudeVLxW, $magnitudeVL, $distance); 
    # Formula to find shortest distance of point with straight line           
    # d (P, L) = |VL x W|/|VL|                                             
    $product = $VL x $W;
    $magnitudeVLxW =  $product->Magnitude();
    $magnitudeVL = $VL->Magnitude();
    $distance = $magnitudeVLxW/$magnitudeVL;
    return $distance;
}

# Description: Finds first and last point on straight line 
# Inputs: PDB file containg points of straight line
# Outputs: First and last point of straight line 2) number of total points on 
#          straight line
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub getBestFitLineCoords    
{                      
    my ($regressionFileHandle) = @_;                                          
    my ($x, $y, $z, @regressionCoords); 
    
    my @regLineCoords = <$regressionFileHandle>;

    foreach my $line (@regLineCoords)
    {
	chomp $line; 
	if ($line =~ /LIN/)
	{
	    my %parsedPDB= parseATOMLine($line);
	    
	    $x = $parsedPDB{x};
	    $y = $parsedPDB{y};
	    $z = $parsedPDB{z};
	    push(@regressionCoords, [$x, $y, $z]);
	}
    }
    return @regressionCoords;                                                 
}  

# Description: Generate figure of epitopes 
# Inputs: PDB file containing co-ordinates of epitope and straight line
# Outputs: A Figure with epitope and straight line
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub makePeptideFigure
{
    my ($PeptidePdbFile, $count) = @_;
    
    # Writing pymol script
    open (my $PYMOL, '>', "Peppymolscript.pml") or die
	"Can not open file"; # opening perl script file 
    my $pdbFileName = basename($PeptidePdbFile, ".pdb");
    
    print {$PYMOL} <<__EOF;
load $PeptidePdbFile
bg_color white
turn x, 90
turn y, 90
turn x, 90
turn y, 90
show cartoon
ray 800,600
png $pdbFileName.png
quit
__EOF
# Sending pymol script to program
`$pymol -c Peppymolscript.pml`;
    unlink "Peppymolscript.pml";
}

# Description: Assigns Secoundary Structure to epitopes from Xmas files
# Inputs: 1) PDB File ID, 2) Antigen chain label, 3) Range of epitope, 
#         4) flag count to number the SS output file
# Outputs: An array containing Secondary structure assignments for epitope
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub SecondaryStrAssignmentToPeptides
{
    my ($pdbFile, $antigenChainLabel, 
	$startRes, $endRes, $count) = @_;
    my ($pdbID, $ext1, @epitopeSS, $ext);
    # Open Output file
    open (my $SS, ">EpitopeSS$count.txt") or
	die "Can not open file $!";
    
    # Xmas files parser to parse secondary structure assignments
#    my $xmastoss = $SFPerlVars::xmastoss;
 #   my $pdbFilePath = "/acrm/data/xmas/pdb/pdb";
 #   $ext = ".xmas";
  #  ($pdbID, $ext1) = split('_',  $pdbFile);
  #  $pdbFile = $pdbFilePath.lc($pdbID).$ext;
  #  $antigenChainLabel = uc($antigenChainLabel);

    my $pdbFilePath = "/acrm/data/pdb/pdb";
    
    $ext = ".ent"; 
    
    ($pdbID, $ext1) = split('_',  $pdbFile);
    $pdbFile = $pdbFilePath.lc($pdbID).$ext;
    $antigenChainLabel = uc($antigenChainLabel);

    
    
    # Parse Secondary structure assignments of only Antigen chain
#   my @antigen =
#	`$xmastoss $pdbFile | grep $antigenChainLabel` ;

    my @antigen = `pdbsecstr $pdbFile | grep $antigenChainLabel`;

      
    my ($insection, $endsection) = 0; 
    # Obtain SS assignments of the given range of epitope
    foreach my $line (@antigen)
    {
	chomp $line;
	# Checks the Antigen chain
	if ($endsection)
	{
	    last;
	}
	if ($line =~ m/$antigenChainLabel$startRes\s+/g)
	{
	    $insection = 1; 
	}
	if ($insection)
	{
	    push (@epitopeSS, $line);
	    print {$SS} $line, "\n"; # Writes on file
	}
	if ($line =~ m/$antigenChainLabel$endRes\s+/g)
	{
	    $endsection=1;
	}
    }
    return @epitopeSS; 
}

# Description: Prints corresponding elements of 2 arrays in one line
# Inputs: 1) 2 Arrays References 2) PDB ID for writing name of PDB 
#         3) file handle 
# Outputs: Write a file combining information from both of the arrays 
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba
sub printArrays
{
    my ($array1Ref, $array2Ref, $pdbID, $ALL) = @_;
    my $count = 1;

    open (my $OUT, ">EpitopeSS_Distance.txt") or die "Can not open";
    my @array1DeRef = @{$array1Ref};
    my @array2DeRef = @{$array2Ref};
    
    my (@eleFromFirst, @eleFromSecond);
    
    print {$ALL} "$pdbID\n";
    print {$OUT} "$pdbID\n";
    
    for( my $rec = 0; $rec < @array1DeRef; $rec++)
    {
	my @eleFromFirst =  @{$array1DeRef[$rec]};
	my @eleFromSecond = @{$array2DeRef[$rec]};
	
	print {$ALL} "Region_$count\n";
	printf {$ALL} "Deviation from beta strand \n", "\n";
	
	print {$OUT} "Region_$count\n";
	printf {$OUT} "Deviation from beta strand \n", "\n";
	
	for (my $i = 0; $i<@eleFromFirst; $i++)
	{
	    printf {$ALL} "%s\t%.3f\n",  
	    $eleFromFirst[$i],$eleFromSecond[$i], "\n\n\n";
	    printf {$OUT} "%s\t%.3f\n",
	    $eleFromFirst[$i],$eleFromSecond[$i], "\n\n";
	}
	print {$ALL} "\n\n";
	$count++;
    }
}

# Description: Finds minimum values in a array
# Inputs: An array with distances of CA with straight line
# Outputs: 1) Position/Index of CA with minimum distance 2) Minimum Distance
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba

sub getLeastDistanceCAlpha
{
    my (@CAalphaDistance) = @_;
    my (%leastCAlpha, $leastDist, $index); 
    
    $leastDist = $CAalphaDistance[0];
    $index = 1; 
    for (my $i = 1; $i <= $#CAalphaDistance; $i++)
    {	    
	if ( $leastDist > $CAalphaDistance[$i])
	{
	    $leastDist =  $CAalphaDistance[$i];
	    $index = $i + 1;
	}
    }
    return ($index, $leastDist);   
}

# Description: Maps Reference point on straight line with constant spacing of
#              3.5 A
# Inputs: 1) 2 vectors; one for straight line and second for point (CA) and
#         first point of straight line 2) closest point index 3) closest point
#         distance 4) PDB file with co-ordinates of point 4) Distances of all
#         points with straight line
# Outputs: A hash with all reference points mapped on line
# Subroutine call/Testing: 
# Date: 03 Nov 2014                                 
# Author: Saba
sub getReferencepoints
{
    my ($coordsCARef, $CAIndex, $closestPointDistance, $P0, $P1, 
	$regressionOF, $ss, $SUMRY, $SSINFO) = @_;
    
    my @coordsCA = @{$coordsCARef};
    my ($refPoint, %refPoints, $move, $numofCA, $forward, $backward, 
	$lineSpacing, $refPointForward, $refPointBackward, $refNumsF, 
	$refNumsB);
    # Getting co-ordinates of closest point 
    my @closestCApoint = @{$coordsCA[$CAIndex-1]}; 
    # Creating Vector for closest point 
    my $P = getVecPointByArrayofPoints (@closestCApoint);

    # Obtaining Direction
    my $VL = $P1 - $P0; # Regression line: P0 is start P1 is end point
    my $W = $P - $P0; # Vector of start point and closest point: 
    # P is closest point

    my $magnitudeVL = $VL->Magnitude();
    my $magnitudeW = $W->Magnitude();

    print {$SUMRY} "Magnitude of best fit line is: $magnitudeVL\n";
    print {$SUMRY} "Magnitude of vector to closest point is: $magnitudeW\n";

    my $lineSegment = getSide ($magnitudeW, $closestPointDistance);
    print {$SUMRY} "The Distance on line from start of line to the point ".
	"(closest to midpoint): $lineSegment\n";
    $move = "CA";
    # Getting First Reference point, point from closest CA
    my $refPointCA = getPointfromDistance($VL, $P0, $lineSegment, 
					  $move);
    print {$SUMRY} "My first reference point is: $refPointCA\n";
    $refPoints{"R$CAIndex"} = $refPointCA; 
  
    # count number of moves in for back tracing and forward tracing
    $numofCA = scalar @coordsCA;
    $forward = $numofCA - $CAIndex; # number of points after first RF
    $backward = $CAIndex - 1; # number of points before first RF
    $refNumsF = $CAIndex + 1; # To define next forward point to the first RF
    $refNumsB = $CAIndex - 1; # To define previous backward point to the first
                              # reference point
   
    # Get percentage of Alpha Helix elements that would be used to decide the
    # spacing of points on the regression line
    my ($alphaPercentage, $betaPercentage, $coilPercentage)
        = getPercentAlpha(@{$ss}); 
    print {$SUMRY} "The percentage of Alpha helices is: $alphaPercentage\n";
    
    if ($alphaPercentage > 60)
    {
	$ss = "alpha"; 
	$lineSpacing = 1.5; # 5.4/3.6 -- alphaSpacing = 1.5
#	print "Alpha Spacing is used\n\n\n";
        print {$SSINFO} "HELIX\n";
    }
    else
    { 
	$ss = "beta";
	$lineSpacing = 3.5; # 7/2 -- beta spacig = 3.5
	print "Beta Spacing is used\n\n\n";
        if ($betaPercentage > 60) {
            print {$SSINFO} "SHEET\n";
        }
        else {
            print {$SSINFO} "COIL\n";    
        }
    }
#######################
    if ($forward)
    {
	$refPoint = $refPointCA; 
	$move = "forward"; 
	for (my $i = 0; $i< $forward; $i++)
	{	    
	    $refPointForward =
		getPointfromDistance ($VL, $refPoint, $lineSpacing, $move);
	    $refPoint = $refPointForward ; 
	    $refPoints{"R$refNumsF"} = $refPoint;  
	    $refNumsF++; 
	}
    }
    if ($backward)
    {
	$refPoint = $refPointCA;
	$move = "backward"; 
	for (my $j = 0; $j< $backward; $j++)
	{
	    $refPointBackward =
		getPointfromDistance ($VL, $refPoint, $lineSpacing, $move);
	    $refPoint = $refPointBackward ;
	    $refPoints{"R$refNumsB"} = $refPoint;
            $refNumsB--;
	}
    }
    # To write the reference points as coorinates/ATOM records of PDB file 
    printAtomeLine ("R", "X", "PON", $regressionOF, %refPoints);
    print {$SUMRY} "The reference points has been written in regression PDB ".
	"as PDB atom lines\n";
#*    makePeptideFigure ($regressionOF, 2);
    return ($ss, %refPoints); 
}
# use pythagorus Formula to get length of side
# Description: 
# Inputs: 
# Outputs: 
# Subroutine call/Testing: 
# Date: 2014                                 
# Author: Saba

sub getSide
{
    my ($hyptenous, $verticalSide) = @_;
    my $base;
    # c^2 (Hypotenous) = a^2 (Base) + b^2 (Vertical Side); 
    # a^2 = c^2 - b^2
    my $product = ($hyptenous ** 2 - $verticalSide ** 2);
    if ($product < 0)
    {
	$product = 0; 
    }
    $base = sqrt ($product);
    return $base; 
}
# Description: 
# Inputs: 
# Outputs: 
# Subroutine call/Testing: 
# Date: 2014                                 
# Author: Saba
sub getPointfromDistance 
{
    my ($VL, $startPoint, $lineSegment, $move) = @_;
    my $magnitudeVL = $VL->Magnitude();
    my $U = $VL/$magnitudeVL;
    my $du = $lineSegment * $U;
    my $refPoint; 
       
    if (($move eq "forward") or ($move eq "CA") or ($move eq "") )
    {
	$refPoint = $startPoint + $du; # Add $du to start point for forward 
    }
    elsif (($move eq "backward"))
    {
	$refPoint = $startPoint - $du;# Substract $du from start point for
    }                                 # backward
    return $refPoint; 
}

# Description: 
# Inputs: 
# Outputs: 
# Subroutine call/Testing: 
# Date: 2014                                 
# Author: Saba
sub getEuclideanDistance
{
    my ($RegionCARef, %refPoints) = @_;
    my $count = 1;
    my @RegionCA = @{$RegionCARef};
    my (@RefPointToCADistance, $distance, $regionCApoint);
    for (my $i = 0; $i < scalar @RegionCA; $i++)
    {
        $regionCApoint = join (",", @{$RegionCA[$i]}); 
	my $key = "R".($i+1);
	$distance = absoluteDistancePoints ($refPoints{$key}, $regionCApoint); 
	push (@RefPointToCADistance, $distance);
    }
    return @RefPointToCADistance; 
}
# Description: 
# Inputs: 
# Outputs: 
# Subroutine call/Testing: 
# Date: 2014                                 
# Author: Saba

sub averageDeviation
{
    my (@RefPointToCADistance) = @_;
    # To add all elements of an array
    my $sum = reduce { $a + $b } @RefPointToCADistance;
    my $average = $sum/((scalar @RefPointToCADistance)); 
    print "Average Deviation :", $average, "\n";
    return $average;
}

# To print co-ordiantes of reference point in PDB format
############################ Distance #####                                   
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                   
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba
sub printAtomeLine
{
    my ($PointKey, $chainName, $resName, $regressionOF, %refPoints) = @_;
    my @ordered_fields = ( 'type', 'serial', 'name', 'altLoc', 'resName',
                           'chainID', 'resSeq','iCode', 'x', 'y', 'z',
                           'occupancy', 'tempFactor', 'element', 'charge');
    my %valueHash = (type => "ATOM", serial => "", name => $chainName, 
		     altLoc => "", resName => $resName, chainID => $chainName,
		     resSeq=>"", iCode => "", x => 1, y => 1, z => 0, 
                     occupancy => 1.00, tempFactor => 1.00, element => "",
		     charge => "");

    open (my $IN, '>>',$regressionOF) or 
	die ("CAN NOT OPEN\n");
    my $count = 1;
    my @atomLines = ();
    foreach my $key (keys %refPoints)
    {
	$valueHash{x} = $refPoints{"$PointKey".$count}->[0];
	$valueHash{y} = $refPoints{"$PointKey".$count}->[1];                   
	$valueHash{z} = $refPoints{"$PointKey".$count}->[2];                  
	$valueHash{serial} = $count;
	$valueHash{resSeq} = $count; 
	$count++;
	my @ordered_values = map {$valueHash{$_}} @ordered_fields;
	my $string = sprintf(  '%-6.6s%5.5s' . ' ' .
			       '%s%1.1s%3.3s %1.1s%4.4s%1.1s      '       
                    .'%8.3f' x 3 .  '%6.2f' x 2 . ' ' x 10 . '%2.2s'         
			       . '%2.2s' ,                                    
                               @ordered_values );   	
	push (@atomLines, $string);
    }
    print {$IN} join("\n" ,@atomLines), "\n\n";
}

# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub absoluteDistancePoints
{
    my ($P, $C) = @_;
    my ($P_x, $P_y, $P_z) = split(/,/, $P);
    my ($C_x, $C_y, $C_z) = split(/,/,$C);
    
    # Formula for absolute distance between 2 points is                   
    # d = sqrt (x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2  
    my $dist = sqrt (($P_x - $C_x) ** 2 +
		     ($P_y - $C_y) ** 2+
		     ($P_z - $C_z) ** 2);
    return $dist; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba
sub reverseArr
{
    my @array;
    push @array, pop @_ while @_;
    return @array; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba
sub getMidPoints
{ 
    # Inputs are vector points
    my ($firstPoint, $lastPoint) = @_;
    my $midPoint = NewVec ( ($firstPoint->[0] + $lastPoint->[0])/2, 
			    ($firstPoint->[1] + $lastPoint->[1])/2,
			    ($firstPoint->[2] + $lastPoint->[2])/2 );
    return $midPoint; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getNearestCAToMidPoint
{
    # Mid point of line
    # CAs of epitope
    my ($midPoint, @coordsCA) = @_;
    my (@MidpointCAdistances, $distance) ; 
    for ( my $i = 0; $i<= $#coordsCA; $i++) 
    {
	$coordsCA[$i] = getVecPointByArrayofPoints (@{$coordsCA[$i]}); 
	$distance = absoluteDistancePoints ($midPoint, $coordsCA[$i]);
	push (@MidpointCAdistances, $distance);
    }
    my ($CAIndex, $closestCADistance) =
	getLeastDistanceCAlpha(@MidpointCAdistances);
    return ($CAIndex, $closestCADistance, @MidpointCAdistances); 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub parseATOMLine {
    my ($ATOM_line) = @_;
    $ATOM_line = pack ( "A81", $ATOM_line );
    my %record
        = ( ATOM => rm_trail( substr($ATOM_line, 0, 6) ),
            serial =>  rm_trail( substr($ATOM_line, 6, 5) ),
            name => rm_trail( substr($ATOM_line, 12, 4) ),
            altLoc => rm_trail( substr($ATOM_line, 16, 1) ),
            resName => rm_trail( substr($ATOM_line, 17, 3) ),
            chainID => rm_trail( substr($ATOM_line, 21, 1) ),
            resSeq => rm_trail( substr($ATOM_line, 22, 4) ),
            iCode => rm_trail( substr($ATOM_line, 26, 1) ),
            x => rm_trail( substr($ATOM_line, 30, 8) ),
	    y => rm_trail( substr($ATOM_line, 38, 8) ),
	    z => rm_trail( substr($ATOM_line, 46, 8) ),                       
            occupancy => rm_trail( substr($ATOM_line, 54, 6) ),             
            tempFactor => rm_trail( substr($ATOM_line, 60, 6) ),
            element => rm_trail( substr($ATOM_line, 76, 2 ) ),
            charge => rm_trail( substr($ATOM_line, 78, 2) ),
	);
    return %record;
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                   
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub rm_trail {
    my $str = shift;
    $str =~ s{\A \s+|\s+ \z}{}gxms;
    return $str;
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getVecPointByArrayofPoints
{
    my (@pointsArray) = @_; 
    my ($P_x, $P_y, $P_z) = @pointsArray; 
    my $point = NewVec($P_x, $P_y, $P_z);
    return $point; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getVectors
{
    my ($startRef, $endRef) = @_; 
    my ($P0, $P1);
    $P0 = getVecPointByArrayofPoints (@{$startRef});                          
    $P1 = getVecPointByArrayofPoints (@{$endRef});  
    my $Vec = $P1 - $P0; 
    return $Vec; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getVectorAngle
{
    my ($VR, $VP) = @_; 
    my $dotProduct = $VR * $VP;
    my ($magVR, $magVP, $cosTheta, $angle, $angleInDegrees); 
    $magVR = $VR->Magnitude();
    $magVP = $VP->Magnitude();
    
    # To find the angle between 2 vectors; 
    # used this formula -- cos (theta) = VR . VP/|VR|*|VP|
    $cosTheta = $dotProduct/($magVR*$magVP); 
    # Inverse of cos (theta) by using acos functiom
    $angle = acos($cosTheta);   
    # convert radians into degree
    $angleInDegrees = rad2deg($angle);
    return $angleInDegrees; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub reverseRegLineDir
{
    my ($regressionFileHandle, $outputFile) = @_;
    my @regPDBFile = <$regressionFileHandle>;
    my (@regArr, @pepArr, @reverseRegArr);
    foreach my $regLine (@regPDBFile)
    {
	if ($regLine =~ /LIN/)
	{
	    push (@regArr, $regLine); 
	}
	else
	{
	    push (@pepArr, $regLine);
	}
    }
    @reverseRegArr =  getRegressionCoordsReversed(@regArr);
    open (my $OUT, '>', $outputFile) or 
	die "$! not open for writing";
    
    print {$OUT} join ("\n", @reverseRegArr), "\n"; 
    print {$OUT} @pepArr; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getResSeq
{
    my ($regressionFileHandle) = @_; 
    my @regPDBFile = <$regressionFileHandle>;
    my (@regArr, @pepArr);
    my (@resSeq, $resSeq); 
    my (@icode, $icode); 
    foreach my $regLine (@regPDBFile)
    {
	chomp $regLine; 
	if ($regLine =~ /LIN/)
	{
	    push (@regArr, $regLine);
	}
        else
        {
            push (@pepArr, $regLine);
        }
    }
    foreach my $atomLine (@pepArr)
    {
	chomp $atomLine;
 	if (isAtomLine($atomLine))
	{
	    my %atomLine = parseATOMLine($atomLine);
	    my $resSeq = $atomLine{resSeq};
	    # To get Insersion code with resSeq
	    $icode = $atomLine{iCode}; 
	    push (@resSeq, $resSeq.$icode);
	}
	else
	{
	    next; 
	}
    }
    # To get just one resSeq out of so many in each ATOM record
    my %tmpHash = map {$_ => 1} @resSeq;
    my @uniqueresSeq = keys %tmpHash; 
    @uniqueresSeq = sort { $a <=> $b } @uniqueresSeq ;
    return @uniqueresSeq; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub isAtomLine
{
    my ($atomLine) = @_;
    if ($atomLine =~ /^(?:ATOM|HETATM)/)
    {
	return 1; 
    }
    else 
    {
	return 0; 
    }
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getPercentAlpha
{
    my (@ss) = @_;
    my $alphaCount = 0; 
    my $betaCount = 0;
    my $coilCount = 0;
    
    foreach my $ssElem (@ss)
    {
	if ($ssElem =~ /(\s+)H|h(\s?)/)
	{
	    $alphaCount++; 
	}
        elsif ($ssElem =~ /(\s+)E|e(\s?)/) {
            $betaCount++;
        }
        elsif ( $ssElem =~ /(\s+)C|c|-(\s?)/) {
            $coilCount++;
        }
    }
    my $alphapercentage = ($alphaCount/scalar @ss) * 100; 
    my $betapercentage = ($betaCount/scalar @ss) * 100;
    my $coilpercentage = ($coilCount/scalar @ss) * 100;
    
    return ($alphapercentage, $betapercentage, $coilpercentage); 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getAllCADistanceWithLine
{
    # I have used first and last reference point on the line of regression
    # to find the vector of straight line 
    my ($startRegRef, $endRegRef, $coordsCARef, $refPointsRef) = @_;     
    my ($P0, $P1, $P); 
    $P0 = getVecPointByArrayofPoints (@{$startRegRef});
    $P1 = getVecPointByArrayofPoints (@{$endRegRef});
    my $RpointVec = $P1-$P0; 
    my @CAdistToLine = calculateDistanceWithLine 
	($P0, $RpointVec, @{$coordsCARef});
    return $RpointVec, @CAdistToLine; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getPointsForCAonLine
{
    my ($RefpointsVec, $RFRef, $coordsCARef, $CAdistToLineRef, 
	$regressionOF) = @_; 
    # $coordsCARef is having arrays of arrays
    open (my $IN, "$$regressionOF") or die "Can not open file";    
    my @coordsCA = @{$coordsCARef}; 
    my @CAdistToLine =  @{$CAdistToLineRef}; 
    # P0 is start of the line
    my $P0 = getVecPointByArrayofPoints (@{$RFRef});
    my ($W, %CApointsOnLine, $CArefPoints, $magnitudeW, $lineSegment); 
    my $count = 0; 
    
    # To find the fall of CA on regression line and map the points for CA on it
    for (my $i = 0; $i <= scalar @{$coordsCARef}-1; $i++)
    {
	my $CA = getVecPointByArrayofPoints (@{$coordsCA[$i]});
	# Vector formed between CA and first point - Hypotenous
	$W = $CA - $P0; 
	$magnitudeW = $W->Magnitude();	
	$lineSegment = getSide($magnitudeW, $CAdistToLine[$i]);
	$count++;
	$CArefPoints = getPointfromDistance($RefpointsVec, $P0,
					    $lineSegment, "");
	$CApointsOnLine{"C".$count} =  $CArefPoints; 
    }
    # To write the co-ordinates of CA points on line as PDB atom lines
    printAtomeLine ("C", "Y", "CAP", $$regressionOF, %CApointsOnLine);
    return %CApointsOnLine; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getPointToPointDistanceOnLine
{
    my ($refPointsOnLineRef, $CApointsOnLineRef) = @_; 
    my %refPointsOnLine = %{$refPointsOnLineRef};
    my %CApointsOnLine = %{$CApointsOnLineRef}; 
    my (@distances, $distance); 
    my $size = keys %refPointsOnLine;
    for (my $i = 1; $i<=$size; $i++)
    {
	$distance = 
	    absoluteDistancePoints ( $refPointsOnLine{"R".$i}, 
				     $CApointsOnLine{"C".$i} );
	print "R".$i." - "."C".$i." = ", $distance, "\n";
	push (@distances, $distance); 	
    }
    return @distances; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getRegressionCoordsReversed
{
    my (@PDB) = @_;
    my (@newPDB, @resSeq, @x, @y, @z, %atomLine, $x, $y, $z); 
    foreach my $atomLine (@PDB)
    {
	%atomLine = parseATOMLine($atomLine);
	my $resSeq = $atomLine{resSeq};
	my $x = $atomLine{x};
	my $y = $atomLine{y};
	my $z = $atomLine{z};
	push (@resSeq, $resSeq);
	push (@x,  $x);
	push (@y,  $y);
	push (@z,  $z);   
    }
    my @ordered_fields = ( 'type', 'serial', 'name', 'altLoc', 'resName',
                           'chainID', 'resSeq','iCode', 'x', 'y', 'z',
		       'occupancy', 'tempFactor', 'element', 'charge');

    my %valueHash = (type => "ATOM", serial => "", name => " O  ",  altLoc 
		     => "",
                     resName => "LIN", chainID => "X", resSeq=>"",
                     iCode => "", x => 1, y => 1, z => 0, occupancy => 1.00,
                     tempFactor => 1.00, element => "", charge => "");

    my $cnt = scalar @resSeq-1; 
    # We need to reverse the co-ordinates leaving the resSeq amd serial num
    # same
    foreach my $rs (@resSeq)
    {
	$valueHash{resSeq} = $rs; 
	$valueHash{serial} = $rs;
	$valueHash{x} = $x[$cnt]; 
	$valueHash{y} = $y[$cnt];
	$valueHash{z} = $z[$cnt];
	my @ordered_values = map {$valueHash{$_}} @ordered_fields;
        my $string = sprintf(  '%-6.6s%5.5s' . ' ' .
                               '%s%1.1s%3.3s %1.1s%4.4s%1.1s   '
                    .'%8.3f' x 3 .  '%6.2f' x 2 . ' ' x 10 . '%2.2s'
                    . '%2.2s' ,                                                
			       @ordered_values );  
	push (@newPDB, $string); 
	$cnt--;  
    }
    return @newPDB; 
}
# Description: The purpose for this subroutine is to find the second closest
# point with the mid point if the last point is the closet point then it 
# shoul look for second closest CA point that should be somewhere in the middle
# Inputs:                                                                    
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getSecondclosestpoint
{
    my ($CAIndex, @MidpointCAdistances) = @_;
    my $sizePep = scalar @MidpointCAdistances; 
    if ($CAIndex == $sizePep)
    {
	splice @MidpointCAdistances, $#MidpointCAdistances, 1;                
    }
    elsif ($CAIndex == 1)
    {
	splice @MidpointCAdistances, 0, 1;
    }
    my ($CAIndex, $closestCADistance) =                                       
        getLeastDistanceCAlpha(@MidpointCAdistances); 
    return ($CAIndex, $closestCADistance); 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub SumOfDistance
{
    my ($beta, @RefPointToCADistance) = @_;
    my ($sum, $mean, $delta);
    # To add all elements of an array      
    $sum = reduce { $a + $b } @RefPointToCADistance;
    print "Summation of Distances is :", $sum, "\n";
    $mean = $sum/((scalar @RefPointToCADistance) ** 2);
    print "Average Distance :", $mean, "\n";
    $delta = $mean - $beta;
    printf "Deviation from Extended beta is : %.3f\n", $delta, "\n";
    return $delta;
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub checkDeviationOnEnds
{
    my (@deviations) = @_;
    my ($firstDev, $lastDev, $avgDev, $lastIndex);
    $lastIndex = $#deviations; 
    $firstDev = $deviations[0];
    $lastDev = $deviations[$#deviations];
    
    # To splice the deviation of last/first 3 residues in case of HOOK
    # and then taking average of the rest of deviations
    if ($firstDev > $lastDev)
    {  
	splice @deviations, 0, 1;
	$avgDev = averageDeviation(@deviations);
	if ($avgDev > 0.5)
	{
	    for ( my $i = 1; $i<3; $i++)
	    {
		splice @deviations, $i-1, 1;
		$avgDev = averageDeviation(@deviations);
		if ($avgDev < 0.5)
		{
		    last; 
		}
	    }
	}
	else 
	{
	    $avgDev = $avgDev;
	}
    }
    else
    {
	splice @deviations, $lastIndex, 1;
	$avgDev = averageDeviation(@deviations);
	my $lastIndex = $#deviations;
	if ($avgDev > 0.5)
        {
	    for ( my $i = 1; $i<3; $i++)
	    {
		my $lastIndex = $#deviations;
		splice @deviations, $lastIndex, 1;
		$avgDev = averageDeviation(@deviations);
		if ($avgDev < 0.5)
                {
                    last;
                }
	    }
        }
        else
	{
            $avgDev = $avgDev;
        }
    }    
    $avgDev = sprintf "%.2f", $avgDev; 
    return ($avgDev);
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub checkAngle
{
    my ($angle) = @_;
    my $shape;  
    if ($angle < 45 )
    {
	$shape = "curved";
    }
    elsif ($angle > 45)
    {
	$shape = "folded";
    }
    return $shape; 
}
# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getContacts
{
    my ($chain, $pdbFile, $sizePep, $ss, $SUMRY, @resSeq) = @_;
    my $length = $sizePep; 
    my ($distance,$maxcontactsLocal); 
    my $ncontacts = 0; 
    my $maxcontactsDistant = 0; 
    my $maxcontactsLocal = 0;
    
    # To check every residue
    my ($countLocal, $countDistant) = 0; 
    my $countL = 0; 
    my $countD = 0;
    for (my $n = 0; $n < $length; $n++)
    {
	# The spacing to check the contacts
	for (my $i = 3; $i < $length - 3; $i++)
	{
	    $ncontacts = 0; 
	    # To count the number of contacts
	    for (my $d = 0; $d < $length; $d++)
	    {
		my $res1 = $n - $d;
		my $res2 = $n + $i + $d;
		if ( ($res1 >= 0) and ($res2 < $length) )
		{
		    $distance =
			getResDist ($resSeq[$res1], $resSeq[$res2], 
				    $chain, $pdbFile);
		    if ($distance)
		    {
			if ($distance <= 4.0)
			{
			    print {$SUMRY} "$n\t$i\t$d\n";
			    print {$SUMRY} "$resSeq[$res1], $resSeq[$res2]\n";
#			    print "$n\t$i\t$d\n";
#                            print "$resSeq[$res1], $resSeq[$res2]\n";
			    print {$SUMRY} "$distance\n";
			    $ncontacts++; 
			}
		    }
		    else
		    {
			next;
		    }
		}
	    }
             # If the length of peptide is more than 12 then it should start 
             # from the spacing of 5 other wise half of length
	    my $contactThershold;  
	    if ($sizePep <= 12)
	    {
		$contactThershold = $sizePep/2; 
	    }
	    else
	    {
		$contactThershold = 5; 
	    }
	    if ($i <= $contactThershold)
	    {
		if ($ncontacts > $maxcontactsLocal)
		{ 
		    $maxcontactsLocal = $ncontacts; 
		}
	    }  
	   elsif ($i > $contactThershold)
	   {
	       if ($ncontacts > $maxcontactsDistant)
	       {
		   $maxcontactsDistant = $ncontacts;
	       }
	   }
	}
#	last; 
    }
    print "Local MaxCOntacts; ", $maxcontactsLocal, "\n";
    print "Distant MaxCOntacts; ", $maxcontactsDistant, "\n";

    print {$SUMRY} "Local MaxCOntacts; ", $maxcontactsLocal, "\n";
    print {$SUMRY} "Distant MaxCOntacts; ", $maxcontactsDistant, "\n";
    my $totalContacts = $maxcontactsLocal + $maxcontactsDistant; 
 
    return ($maxcontactsLocal, $maxcontactsDistant, $totalContacts)
}

# Description:                                                                
# Inputs:                                                                     
# Outputs:                                                                    
# Subroutine call/Testing:                                                    
# Date: 2014                                                                  
# Author: Saba 
sub getResDist
{
    my ($res1, $res2, $chain, $pdbFile) = @_;
    my $resdist = $SFPerlVars::resdist;
    $res1 = $chain.".".$res1;
    $res2 = $chain.".".$res2;
    my $distance = `$resdist $res1 $res2 ../$pdbFile`;
    return $distance; 
}
1;
