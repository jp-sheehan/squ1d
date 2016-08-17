#!/bin/bash

# defaults
program="gnuplot"
dir="./"

nspec=$(grep -n "SPECIES" SolverInput.inp | cut -c1)
linenum=$((${nspec}+3))
ionspec=$(sed "${linenum}q;d" SolverInput.inp | cut -d" " -f10)

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

pusher=$(fgrep "FIELD_LOCATION" SolverType.inp | head -1 | grep -o "1" | tr -d '\n')
if [[ "$pusher" == "11" ]]; then
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
