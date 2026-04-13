#!/bin/bash

export path_xyz="../../Data/MD_Data/AgCo/XYZ"
export electron_energy=200
export nz=20

export process_id=$1

mkdir -p tmp_$process_id/xyz

$ovitos cluster_separation.py   # creates xyz files in ./xyz as {nanoparticle_id}_{id_sim}.xyz

$ovitos particle_analyze.py > /dev/null 2>&1    # creates data.dat in the current directory

./create_msa_prm.py

for xyz_file in tmp_$process_id/xyz/*.xyz
do
    ./mk_image.py $xyz_file
done
