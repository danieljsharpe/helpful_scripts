#!/bin/bash

# bash script to process an "Epath" format file, and use PATHSAMPLE to extract minima and transition states
# along the path and write a .crd (AMBER) format trajectory file, that can be read by e.g. VMD
# make sure you have a PATHSAMPLE executable named "PATHSAMPLE" in your path!
# pathdata needs to initially have both "EXTRACTMIN" and "EXTRACTTS" lines
# also need the Fortran executable "threecol_2_crd" for converting between PATHSAMPLE output and AMBER .crd formats

# Daniel J. Sharpe, Mar 2019

# set names of input and output files, and number of atoms!
Epath_fname=Epath
natoms=334
crd_fname=Epath.mdcrd

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
        rm extractedmin.rst
        # PATHSAMPLE dumps coordinates in plain 3-column x-y-z format, change to AMBER .crd format
        $script_dir/threecol_2_crd $natoms >> $crd_fname

        cp pathdata pathdata.tmp
        sed 's/! EXTRACTTS/EXTRACTTS/' pathdata.tmp > pathdata.tmp2
        sed 's/EXTRACTMIN/! EXTRACTMIN/' pathdata.tmp2 > pathdata
    elif [ "$ts" = true ]
    then
        echo "    sp is a ts"
        min=true
        ts=false
        cp pathdata pathdata.tmp
        sed 's/.*EXTRACTTS.*/EXTRACTTS '$sp_id'/' pathdata.tmp > pathdata.tmp2
        sed 's/.*EXTRACTMIN.*/! EXTRACTMIN/' pathdata.tmp2 > pathdata

        PATHSAMPLE > PS_extract.out
        mv extractedts extractedmin
        rm extractedts.rst
        # PATHSAMPLE dumps coordinates in plain 3-column x-y-z format, change to AMBER .crd format
        $script_dir/threecol_2_crd $natoms >> $crd_fname

        cp pathdata pathdata.tmp
        sed 's/! EXTRACTMIN/EXTRACTMIN/' pathdata.tmp > pathdata.tmp2
        sed 's/EXTRACTTS/! EXTRACTTS/' pathdata.tmp2 > pathdata
    fi
    rm extractedmin
done < $Epath_fname
