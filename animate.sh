#!/bin/bash

# defaults
program="gnuplot"
dir="./"
ionspec="ARGON"

# inputs

while getopts ":d:i:p:" opt; do
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
   script="~/squ1d/scripts/animateplasma.gp"
fi

for p in "${parameters[@]}"; do
   #script="./plot/${program}_${p}.gp"
   vars="maxnum=$number; param='$p'; loc='$loc'; dir='$dir'; ionspec='$ionspec'"
   cmd="$program -e \"$vars\" $script"
   eval $cmd #run script in this terminal
   # gnome-terminal -e sh $script #run script from new terminal (to have multiple windows)
   gifview "$dir$p.gif"
done
