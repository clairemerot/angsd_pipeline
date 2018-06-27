#!/bin/bash

pop_file="./02_info/pop.txt"
num_pops=$(wc -l "$pop_file" | cut -d " " -f 1)

for i in $(seq $num_pops)
do
	pop1=$(cat "$pop_file" | head -"$i" | tail -1)
	for j in $(seq $[ $i + 1 ] $num_pops)
	do
		pop2=$(cat "$pop_file" | head -"$j" | tail -1)
		echo $pop1 $pop2
		
		# Do your stuff
	done
done
