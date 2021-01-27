#!/bin/bash

# bash script to process an "Epath" format file, and use PATHSAMPLE to extract minima and transition states
# along the path and write a .xyz format trajectory file, that can be read by e.g. VMD
# make sure you have a PATHSAMPLE executable named "PATHSAMPLE" in your path!
# pathdata needs to initially have both "EXTRACTMIN" and "EXTRACTTS" lines
# this version of the script is for processing structures from HiRE, ID'ing beads according to base ID
# need the Python script "plain_2_xyz_hire.py" for converting between PATHSAMPLE (plain) output and .xyz formats

# Daniel J. Sharpe, Mar 2019

# set names of input and output files, and number of atoms!
Epath_fname=Epath_3_U
xyz_fname=Epath_3_U.xyz

script_dir=`dirname "$0"`

rm $crd_fname
echo "default_name" >> $crd_fname
echo "I am going to loop over the Epath file named: $Epath_fname"

min=true
ts=false
while read line; do
    sp_id=$(echo $line | awk '{print $3}')
    echo "sp_id: $sp_id"
    if [ "$min" = true ]
    then
        echo "    sp is a minimum"
        min=false
        ts=true
        cp pathdata pathdata.tmp
        sed 's/.*EXTRACTMIN.*/EXTRACTMIN '$sp_id'/' pathdata.tmp > pathdata.tmp2
        sed 's/.*EXTRACTTS.*/! EXTRACTTS/' pathdata.tmp2 > pathdata

        PATHSAMPLE > PS_extract.out
        # PATHSAMPLE dumps coordinates in plain 3-column x-y-z format, change to .xyz format
        python $script_dir/plain_2_xyz_hire.py extractedmin
        cat extractedmin.xyz >> $xyz_fname

        cp pathdata pathdata.tmp
        sed 's/! EXTRACTTS/EXTRACTTS/' pathdata.tmp > pathdata.tmp2
        sed 's/EXTRACTMIN/! EXTRACTMIN/' pathdata.tmp2 > pathdata
        rm extractedmin extractedmin.xyz
    elif [ "$ts" = true ]
    then
        echo "    sp is a ts"
        min=true
        ts=false
        cp pathdata pathdata.tmp
        sed 's/.*EXTRACTTS.*/EXTRACTTS '$sp_id'/' pathdata.tmp > pathdata.tmp2
        sed 's/.*EXTRACTMIN.*/! EXTRACTMIN/' pathdata.tmp2 > pathdata

        PATHSAMPLE > PS_extract.out
        # PATHSAMPLE dumps coordinates in plain 3-column x-y-z format, change to AMBER .crd format
        python $script_dir/plain_2_xyz_hire.py extractedts
        cat extractedts.xyz >> $xyz_fname

        cp pathdata pathdata.tmp
        sed 's/! EXTRACTMIN/EXTRACTMIN/' pathdata.tmp > pathdata.tmp2
        sed 's/EXTRACTTS/! EXTRACTTS/' pathdata.tmp2 > pathdata
        rm extractedts extractedts.xyz
    fi
done < $Epath_fname
