#!/bin/sh
# usage: svn st | x 2 | xargs rm
# usage 1: sh get_col2.sh 2 infile > outfile # to print col 2 of infile to outfile
# usage 2: cat infile | sh get_col2.sh 2 > outfile # does same as above
col=$1
shift
awk -v col="$col" '{print $col}' "${@--}"
