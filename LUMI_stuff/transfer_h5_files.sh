#!/bin/bash

parent_dir="/scratch/project_465001379/varverak/hdf5/front_initial_eta16_propagated"
target_dir="/PIC_Plasma"

for file in "$parent_dir"/*; do
	filename="$(basename -- $file)"
	dbxcli put "$file" "$target_dir/$filename"
	#echo "$file will be copied to $target_dir/$filename"
done
