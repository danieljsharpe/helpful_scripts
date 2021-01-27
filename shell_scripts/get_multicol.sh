#!/bin/bash

awk '{ print $4, $5 }' infile >> outfile
