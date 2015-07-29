#!/bin/bash                                                                    
#*************************************************************************      
#                                                                               
#   Program:    processAntibodyMartin
#   File:       processAntibodyMartin.sh
#                                                                               
#   Version:    V1.1                                                            
#   Date:       22.04.14                                                        
#   Function:   This is an automatic pipeline to process the Data for AbDb and 
#               It generates data just for Martin numbering scheme 
#   Usage:      ./processAntibodyMartin.sh
#                                                                               
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
#*************************************************************************

function makeNRdata
{
mkdir -p "NR_"$proAntigen
cd $proAntigen  # runs perl program for non-redundant data   
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$proAntigen".txt" "NR_"$proAntigen
cd ..

mkdir -p "NR_"$free
cd $free
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$free".txt" "NR_"$free
cd ..

mkdir -p "NR_"$npAntigen
cd $npAntigen
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$npAntigen".txt" "NR_"$npAntigen
cd ..

mkdir -p "NR_"$Combined
cd $Combined
perl ~/scripts/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$Combined.txt "NR_"$Combined
cd ..

}

function process
{
 # parameter
free=$1
proAntigen=$2
npAntigen=$3
Combined=$4

Redundant="Redundant_files"
mkdir -p $Combined
mkdir -p $Redundant
# copy all the complexes in Combined directory
cp ./$free/* ./$Combined 
cp ./$proAntigen/* ./$Combined
cp ./$npAntigen/* ./$Combined

cd $free
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

cd $npAntigen
perl ~/scripts/bin/getRedundantAntibodyClusters.pl; 
mv *.txt ../$Redundant
cd ..

cd $proAntigen
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

cd $Combined
perl ~/scripts/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

makeNRdata
data="Data"

mkdir -p $data
mv $free ./$data
mv $proAntigen ./$data
mv $npAntigen ./$data
mv $Combined ./$data
mv "NR_"* ./$data

} # Function ends
##############################

function countFiles
{
    dir=`pwd`
    data="Data"
    
processed_proAntigenAB=`ls $dir/$data/AntibodyAntigen_Martin | wc -l`
processed_antibody=`ls $dir/$data/FreeAntibody_Martin | wc -l`
processed_haptenAB=`ls $dir/$data/AntibodyHapten_Martin | wc -l`
processed_combinedAB=`ls $dir/$data/CombinedAb_Martin | wc -l`
NR_proAntigenAB=`ls $dir/$data/NR_AntibodyAntigen_Martin | wc -l`
NR_Antibody=`ls $dir/$data/NR_FreeAntibody_Martin | wc -l`
NR_HaptenAB=`ls $dir/$data/NR_AntibodyHapten_Martin | wc -l`
NR_CombinedAB=`ls $dir/$data/NR_CombinedAb_Martin | wc -l`

processed_proAntigenLG=`ls $dir/$data/LightAntigen_Martin | wc -l`
processed_light=`ls $dir/$data/LightChain_Martin | wc -l`
processed_haptenLG=`ls $dir/$data/LightHapten_Martin | wc -l`
processed_combinedLG=`ls $dir/$data/CombinedLg_Martin | wc -l`
NR_proAntigenLG=`ls $dir/$data/NR_LightAntigen_Martin | wc -l`
NR_light=`ls $dir/$data/NR_LightChain_Martin | wc -l`
NR_HaptenLG=`ls $dir/$data/NR_LightHapten_Martin | wc -l`
NR_CombinedLG=`ls $dir/$data/NR_CombinedLg_Martin | wc -l`

processed_proAntigenHV=`ls $dir/$data/HeavyAntigen_Martin | wc -l`
processed_heavy=`ls $dir/$data/HeavyChain_Martin | wc -l`
processed_haptenHV=`ls $dir/$data/HeavyHapten_Martin | wc -l`
processed_combinedHV=`ls $dir/$data/CombinedHv_Martin | wc -l`
NR_proAntigenHV=`ls $dir/$data/NR_HeavyAntigen_Martin | wc -l`
NR_heavy=`ls $dir/$data/NR_HeavyChain_Martin | wc -l`
NR_HaptenHV=`ls $dir/$data/NR_HeavyHapten_Martin | wc -l`
NR_CombinedHV=`ls $dir/$data/NR_CombinedHv_Martin | wc -l`

    
echo "$processed_proAntigenAB,$processed_antibody,$processed_haptenAB,$processed_combinedAB,$NR_proAntigenAB,$NR_Antibody,$NR_HaptenAB,$NR_CombinedAB,$processed_proAntigenLG,$processed_light,$processed_haptenLG,$processed_combinedLG,$NR_proAntigenLG,$NR_light,$NR_HaptenLG,$NR_CombinedLG,$processed_proAntigenHV,$processed_heavy,$processed_haptenHV,$processed_combinedHV,$NR_proAntigenHV,$NR_heavy,$NR_HaptenHV,$NR_CombinedHV"
                                                                
}


function Split 
{
pattern=$1;
String=`grep $pattern masterlog-a.log` # e.g: returns complexes=400
parameter=$(echo $String | cut -f1 -d=) # splits the string and returns number
count=$(echo $String | cut -f2 -d=) # after =, the delimiter
}
##############################


function compress
{
    for directory in  *
    do
        dirname=$(basename "$directory")
        tar -jcvf $dirname.tar.bz2 $dirname
    done 
}


# The main program begins here
function runProg
{
scheme=$1
schemeFlag=$2
# Running get_antibody_complex program for 3 numbering schemes
echo "get_antibody_complex program is running for $scheme numbering";
free="FreeAntibody_"$scheme
proAntigen="AntibodyAntigen_"$scheme
npAntigen="AntibodyHapten_"$scheme
Combined="CombinedAb_"$scheme
perl ~/scripts/bin/processAntibodyPDBs.pl $schemeFlag $3
process $free $proAntigen $npAntigen $Combined

free="LightChain_"$scheme
proAntigen="LightAntigen_"$scheme
npAntigen="LightHapten_"$scheme
Combined="CombinedLg_"$scheme
process $free $proAntigen $npAntigen $Combined

free="HeavyChain_"$scheme
proAntigen="HeavyAntigen_"$scheme
npAntigen="HeavyHapten_"$scheme
Combined="CombinedHv_"$scheme
process $free $proAntigen $npAntigen $Combined

}

############################################
echo "Main program running";
scheme="Kabat"
schemeFlag="-k"
runProg $scheme $schemeFlag $1

scheme="Chothia"
schemeFlag="-c"
runProg $scheme $schemeFlag $1

scheme="Martin"
schemeFlag="-a"
runProg $scheme $schemeFlag $1

Result=$(countFiles)

function get_rtrn(){
    echo `echo $1|cut --delimiter=, -f $2`
}


processed_proAntigenAB=`get_rtrn $Result 1`
processed_antibody=`get_rtrn $Result 2`
processed_haptenAB=`get_rtrn $Result 3`
processed_combinedAB=`get_rtrn $Result 4`
NR_proAntigenAB=`get_rtrn $Result 5`
NR_Antibody=`get_rtrn $Result 6`
NR_HaptenAB=`get_rtrn $Result 7`
NR_CombinedAB=`get_rtrn $Result 8`

processed_proAntigenLG=`get_rtrn $Result 9`
processed_light=`get_rtrn $Result 10`
processed_haptenLG=`get_rtrn $Result 11`
processed_combinedLG=`get_rtrn $Result 12`
NR_proAntigenLG=`get_rtrn $Result 13`
NR_light=`get_rtrn $Result 14`
NR_HaptenLG=`get_rtrn $Result 15`
NR_CombinedLG=`get_rtrn $Result 16`

processed_proAntigenHV=`get_rtrn $Result 17`
processed_heavy=`get_rtrn $Result 18`
processed_haptenHV=`get_rtrn $Result 19`
processed_combinedHV=`get_rtrn $Result 20`
NR_proAntigenHV=`get_rtrn $Result 21`
NR_heavy=`get_rtrn $Result 22`
NR_HaptenHV=`get_rtrn $Result 23`
NR_CombinedHV=`get_rtrn $Result 24`

echo $processed_proAntigenAB
echo $processed_antibody
echo $processed_haptenAB
echo $processed_combinedAB
echo $NR_proAntigenAB
echo $NR_Antibody
echo $NR_HaptenAB
echo $NR_CombinedAB

echo $processed_proAntigenLG
echo $processed_light
echo $processed_haptenLG
echo $processed_combinedLG
echo $NR_proAntigenLG
echo $NR_light
echo $NR_HaptenLG
echo $NR_CombinedLG

echo $processed_proAntigenHV
echo $processed_heavy
echo $processed_haptenHV
echo $processed_combinedHV
echo $NR_proAntigenHV
echo $NR_heavy
echo $NR_HaptenHV
echo $NR_CombinedHV


log="Logs"
mkdir -p $log
mv *.log ./$log 
mv *.list ./$log

mv -f $Redundant ./$data
cd $data/$Redundant

# To sort all the files in directory
shopt -s nullglob
filearray=( * ) # Reading directory files into an array
 for i in "${filearray[@]}"
 do
     sort $i -o $i
 done
cd ..

compress

exit;

# The following bit of code finds stats of processed PDBs from masterlog.log 
cd $log
pattern="AntibodyAntigen="
Split $pattern
AntibodyAntigen=$count

pattern="Free_Antibody="
Split $pattern
Freeantibody=$count

pattern="AB-Hapten"
Split $pattern
ABhapten=$count

pattern="CompleteAntibody_Dataset="
Split $pattern
completeAntibodyDataset=$count

pattern="Light-Antigen="
Split $pattern
lightAntigen=$count

pattern="Bence-Jones="
Split $pattern
bensJones=$count

pattern="Light-Hapten="
Split $pattern
LGhapten=$count

pattern="CompleteLight_Dataset="
Split $pattern
completeLightDataset=$count

pattern="Heavy-Antigen="
Split $pattern
heavyAntigen=$count

pattern="Camelids="
Split $pattern
camelids=$count

pattern="Heavy-Hapten="
Split $pattern
HVhapten=$count

pattern="CompleteHeavy_Dataset="
Split $pattern
completeHeavyDataset=$count


pattern="Fc="
Split $pattern
fc=$count

pattern="Failed="
Split $pattern
kabatFailed=$count

pattern="Superseded="
Split $pattern
superseded=$count

pattern="CDR-Error="
Split $pattern
cdrError=$count

cp Kabat_Failed.list ../$data
cp ../header.dat ../$data
cd ../$data

# To merge 2 consective lines into one
awk 'NR%2{printf $0" ";next;}1' header.dat >headerProcessed.dat
cd ..

#cd ..
bash ~/scripts/bin/statsProcessed.sh $complex $processed_complexes $NR_Complex $hapten $processed_hapten $NR_Hapten $antibody $processed_antibody $NR_Antibody $completeDataset $processed_combined $NR_Combined >stats_processed.tt
bash ~/scripts/bin/statsUnprocessed.sh $bensJones $lightAntigen $camelids $heavyAntigen $fc $kabatFailed $superseded >stats_unprocessed.tt

cd $data
mkdir -p NR_Martin_merged
cd NR_Complex_Martin
cp *.pdb ../NR_Martin_merged
cd ..

cd NR_Antibody_Martin
cp *.pdb ../NR_Martin_merged
cd ..

cd NR_Hapten_Martin
cp *.pdb ../NR_Martin_merged
cd ..

# To copy the PDB codes from the directory  
cd NR_Martin_merged
ls | grep "_" | cut -f1 -d. > ../Merged.txt
cd ..

# This shell command finds diffrence between Combined.txt and Merged.txt
comm -13 <(sort Combined.txt) <(sort Merged.txt) >difference.txt
# This script produces list of antibody PDBs that are present in both forms  
# i.e Free or Complexed 
perl ~/scripts/bin/FreeComplexedAntibody.pl

cd ..
exit;

# This shell script does not update the website...
