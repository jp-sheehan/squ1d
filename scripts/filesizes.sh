#!/bin/bash

filenames=($(ls | egrep [do][ua]t$))
output="FILENAME LINES SIZE SIZE/LINE\n"

for f in "${filenames[@]}"; do
   lines=$(wc -l $f | cut -d' ' -f1)
   size=$(ls -alF $f | cut -d' ' -f5)
   if [ $lines -eq "0" ] ; then
      sizeperline=0
   else
      sizeperline=$((size/lines))
   fi
   output="${output}$f $lines $size $sizeperline\n"
done

echo -en $output | column -t
