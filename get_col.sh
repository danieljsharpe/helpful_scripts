#!/bin/bash

file_old=$1
file_new=$2

while read c1 c2; do echo $c2; done < $file_old > $file_new
