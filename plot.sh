#!/bin/bash

# defaults
program="gnuplot"
dir="./"

# inputs
while getopts ":d:p:n:" opt; do
   case $opt in
      d)
         # set the directory in which to search
         dir=$OPTARG
         ;;
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

number=$(ls -t $dir| grep restart | head -1 | grep -o '[0-9]\+') # most recent number
shift $((OPTIND-1))
parameters=("$@")

if [ $program == "gnuplot" ]
then
   script="./plot/plasmaplot.gp"
fi

for p in "${parameters[@]}"; do
   #script="./plot/${program}_${p}.gp"
   vars="fnum=$number; param='$p'; dir='$dir'"
   cmd="$program -e \"$vars\" $script"
   eval $cmd #run script in this terminal
   # gnome-terminal -e sh $script #run script from new terminal (to have multiple windows)
done
