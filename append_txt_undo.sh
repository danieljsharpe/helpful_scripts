#!/bin/bash

n_Epath_files=150

for i in $(seq 1 1 $n_Epath_files)
do
    mv Epath.$i.txt Epath.$i
done
