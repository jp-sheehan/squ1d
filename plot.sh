#!/bin/bash

# defaults
program="gnuplot"

while getopts ":p:" opt; do
   case $opt in
      p)
         # Set which program should be used to plot the data
         program=$OPTARG
         ;;
      \?)
         echo "Invalid option: -$OPTARG" >&2
         ;;
      :)
         echo "Option -$OPTARG requires an argument." >&2
         exit 1
         ;;
   esac
done

shift $((OPTIND-1))
parameters=("$@")

for p in "${parameters[@]}"; do
   script="./plot/${program}_${p}.gp"
   echo $program
   echo $script
   cmd="$program $script"
   eval $cmd #run script in this terminal
   # gnome-terminal -e sh $script #run script from new terminal (to have multiple windows)
done
