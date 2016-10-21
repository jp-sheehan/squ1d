#!/bin/bash

data=$(ls | fgrep .dat)
output=$(ls | fgrep .out)
input=$(ls | fgrep .inp)

directory="$1";
mkdir $directory

mv -t $directory $data
mv -t $directory $output
cp -t $directory $input
