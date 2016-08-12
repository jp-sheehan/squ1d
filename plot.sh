#!/bin/bash

# defaults
program="gnuplot"
dir="./"
ionspec="ARGON"
number=-1

# inputs

while getopts ":d:i:p:n:" opt; do
   case $opt in
      d)
         # set the directory in which to search
         dir=$OPTARG
         lastchar="${dir: -1}"
         if [ "$lastchar" != "/" ]; then
            dir="$dir/"
         fi
         ;;
      i)
         # set the ion specie
         ionspec=$OPTARG
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

if [ "$number" -eq "-1" ]; then
   number=$(ls -t $dir| grep restart | head -1 | grep -o '[0-9]\+') # most recent number
fi

pusher=$(fgrep "PARTICLE_MOVER" ${dir}SolverType.inp | head -1)
if [[ "$pusher" == *"SIMPLE"* ]]; then
   loc="p"
elif [[ "$pusher" == *"BORIS"* ]]; then
   loc="p"
elif [[ "$pusher" == *"Q1D"* ]]; then
   loc="p"
else
   loc="c"
fi

shift $((OPTIND-1))
parameters=("$@")

if [ $program == "gnuplot" ]
then
   script="~/squ1d/scripts/plotplasma.gp"
fi

for p in "${parameters[@]}"; do
   #script="./plot/${program}_${p}.gp"
   vars="fnum=$number; param='$p'; dir='$dir'; loc='$loc'; ionspec='$ionspec'"
   cmd="$program -e \"$vars\" $script"
   eval $cmd #run script in this terminal
   # gnome-terminal -e sh $script #run script from new terminal (to have multiple windows)
done
