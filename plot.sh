#!/bin/bash

# defaults
program="gnuplot"
number=$(ls -t | grep restart | head -1 | grep -o '[0-9]\+') # most recent number

# inputs
while getopts ":p:n:" opt; do
   case $opt in
      p)
         # Set which program should be used to plot the data
         program=$OPTARG
         ;;
      n)
         # data set number
         number=$OPTARG
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

if [ $program == "gnuplot" ]
then
   script="./plot/plasmaplot.gp"
fi

for p in "${parameters[@]}"; do
   #script="./plot/${program}_${p}.gp"
   vars="fnum=$number; param='$p'"
   cmd="$program -e \"$vars\" $script"
   eval $cmd #run script in this terminal
   # gnome-terminal -e sh $script #run script from new terminal (to have multiple windows)
done
