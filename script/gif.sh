#!/usr/bin/env bash

## Collect all png files in the files array
rm *.gif
files=(`ls * | sort -nk1,1`)
#echo ${files[@]}
## How many should be done at once
batch=50
j=0
## Read the array in batches of $batch
for (( i=0; $i<${#files[@]}; i+=$batch ))
do
    ## Convert this batch
    convert -delay 5 -loop 0 `echo "${files[@]:$i:$batch}" | sort -nk1,1` $j.animated.gif &
    ((j++))
done

## Now, merge them into a single file
#convert  `ls *.animated.gif | sort -nk1,1` all.gif