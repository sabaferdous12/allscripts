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
perl ~/allscript/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$proAntigen".txt" "NR_"$proAntigen
cd ..

mkdir -p "NR_"$free
cd $free
perl ~/allscript/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$free".txt" "NR_"$free
cd ..

mkdir -p "NR_"$npAntigen
cd $npAntigen
perl ~/allscript/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$npAntigen".txt" "NR_"$npAntigen
cd ..

mkdir -p "NR_"$Combined
cd $Combined
perl ~/allscript/bin/constructNonRedundantData.pl ../$Redundant/"Redundant_"$Combined.txt "NR_"$Combined

#ls *.pdb | grep "_" | cut -f1 -d. ../>"NR_"$Combined".txt"
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
perl ~/allscript/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

cd $npAntigen
perl ~/allscript/bin/getRedundantAntibodyClusters.pl; 
mv *.txt ../$Redundant
cd ..

cd $proAntigen
perl ~/allscript/bin/getRedundantAntibodyClusters.pl;
mv *.txt ../$Redundant
cd ..

cd $Combined
perl ~/allscript/bin/getRedundantAntibodyClusters.pl;
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
# This function puts together all the data (Complete Antibody, Light and
# Heavy chains) for each numbering scheme
function combineData
{
    scheme=$1;
    Redundant="Redundant_files"

    mkdir "ALL_"$scheme
    cp ./"CombinedAb_"$scheme/* ./ALL_$scheme
    cp ./"CombinedLg_"$scheme/* ./ALL_$scheme
    cp ./"CombinedHv_"$scheme/* ./ALL_$scheme

    cd "ALL_"$scheme
    perl ~/allscript/bin/getRedundantAntibodyClusters.pl;
    mv *.txt ../../$Redundant
    cd ..
}

function countFiles
{
    dir=`pwd`
    data="Data"
    cd $data
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

scheme="Kabat"
combineData $scheme
scheme="Chothia"
combineData $scheme
scheme="Martin"
combineData $scheme

mkdir -p  NR_Merged_Antibody
cp ./NR_AntibodyAntigen_Martin/* ./NR_Merged_Antibody
cp ./NR_AntibodyHapten_Martin/* ./NR_Merged_Antibody
cp ./NR_FreeAntibody_Martin/* ./NR_Merged_Antibody
cd ./NR_Merged_Antibody
`ls *.pdb | grep "_" | cut -f1 -d. >../NR_Merged_Antibody.txt`
cd ..

mkdir -p NR_Merged_Light
cp ./NR_LightAntigen_Martin/* ./NR_Merged_Light
cp ./NR_LightHapten_Martin/* ./NR_Merged_Light
cp ./NR_LightChain_Martin/* ./NR_Merged_Light

cd ./NR_Merged_Light
`ls *.pdb | grep "_" | cut -f1 -d. >../NR_Merged_Light.txt`
cd ..


mkdir -p NR_Merged_Heavy
cp ./NR_HeavyAntigen_Martin/* ./NR_Merged_Heavy
cp ./NR_HeavyHapten_Martin/* ./NR_Merged_Heavy
cp ./NR_HeavyChain_Martin/* ./NR_Merged_Heavy

cd ./NR_Merged_Heavy
`ls *.pdb | grep "_" | cut -f1 -d. >../NR_Merged_Heavy.txt`
cd ..

cd ./NR_CombinedAb_Martin
`ls *.pdb | grep "_" | cut -f1 -d. >../NR_Combined_Antibody.txt`
cd ..

cd ./NR_CombinedLg_Martin
`ls *.pdb | grep "_" | cut -f1 -d. >../NR_Combined_Light.txt`
cd ..

cd ./NR_CombinedHv_Martin
`ls *.pdb | grep "_" | cut -f1 -d. >../NR_Combined_Heavy.txt`
cd ..

# This shell command finds diffrence between Combined.txt and Merged.txt
comm -13 <(sort NR_Combined_Antibody.txt) <(sort NR_Merged_Antibody.txt) >difference_Antibody.txt

comm -13 <(sort NR_Combined_Light.txt) <(sort NR_Merged_Light.txt) >difference_Light.txt

comm -13 <(sort NR_Combined_Heavy.txt) <(sort NR_Merged_Heavy.txt) >difference_Heavy.txt

# Create list of free antibody and complexed antibodies 
perl ~/allscript/bin/FreeComplexedAntibody.pl -a ./difference_Antibody.txt ../Redundant_files/Redundant_CombinedAb_Martin.txt

# Create list of free Bence-Jones and complexed Bence-jones
perl ~/allscript/bin/FreeComplexedAntibody.pl -l ./difference_Light.txt ../Redundant_files/Redundant_CombinedLg_Martin.txt

# Create list of free antibody and complexed antibodies
perl ~/allscript/bin/FreeComplexedAntibody.pl -h ./difference_Heavy.txt ../Redundant_files/Redundant_CombinedHv_Martin.txt

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
    array=(*/)
    for directory in  "${array[@]}"
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
perl ~/allscript/bin/processAntibodyPDBs.pl $schemeFlag $3
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

function get_rtrn(){
        echo `echo $1|cut --delimiter=, -f $2`
        }


############################################
echo "Main program running";
scheme="Kabat"
schemeFlag="-k"
runProg $scheme $schemeFlag $1
exit; 

scheme="Chothia"
schemeFlag="-c"
runProg $scheme $schemeFlag $1

scheme="Martin"
schemeFlag="-a"
runProg $scheme $schemeFlag $1

Result=$(countFiles)

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
mv *.txt ../$log 
compress
cd ..
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

bash ~/allscript/bin/statsProcessed.sh $AntibodyAntigen $processed_proAntigenAB $NR_proAntigenAB $ABhapten $processed_haptenAB $NR_HaptenAB $Freeantibody $processed_antibody $NR_Antibody $completeAntibodyDataset $processed_combinedAB $NR_CombinedAB $lightAntigen $processed_proAntigenLG $NR_proAntigenLG $LGhapten $processed_haptenLG $NR_HaptenLG $bensJones $processed_light $NR_light $completeLightDataset $processed_combinedLG $NR_CombinedLG $heavyAntigen $processed_proAntigenHV $NR_proAntigenHV $HVhapten $processed_haptenHV $NR_HaptenHV $camelids $processed_heavy $NR_heavy $completeHeavyDataset $processed_combinedHV $NR_CombinedHV >stats_processed.tt

bash ~/allscript/bin/statsUnprocessed.sh $fc $kabatFailed $cdrError $superseded >stats_unprocessed.tt

exit

