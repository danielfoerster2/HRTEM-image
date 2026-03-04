#!/bin/bash

export path_xyz="/home/victor.barrere@crmd.cnrs-orleans.fr/Documents/Data/Data_xyz/"
export path_processed="/home/victor.barrere@crmd.cnrs-orleans.fr/Documents/Data/Data_processed_test/"
export path_new_xyz=$path_processed"XYZ/"

mkdir -p $path_new_xyz
rm -f $path_processed/data.dat

$ovitos cluster_separation.py
$ovitos particle_analyze.py

