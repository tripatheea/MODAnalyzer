#!/bin/bash

output_file=$(realpath $2)

> ${output_file}

total=0
serial=0
cd $1
for i in *; do
	x=$i
	filename=${x:(0):(-4)} 
	number=$(grep -o BeginEvent $i | wc -l)
	echo $filename $number
	echo $filename $number >> $output_file
	total=$((total + number))
	serial=$((serial + 1))
done
echo "Total" $serial $total
echo "Total" $serial $total >> $output_file

