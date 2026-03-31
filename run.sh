#!/bin/bash

export path_xyz="../Data/MD_data/AgCo/XYZ"
export electron_energy=300
export nz=10

export process_id=$1

cd tmp_$process_id
mkdir -p ./xyz

$ovitos cluster_separation.py   # creates xyz files in ./xyz as {nanoparticle_id}_{id_sim}.xyz
$ovitos particle_analyze.py     # creates data.dat in the current directory

./create_msa_prm.py

for xyz_file in ./xyz/*.xyz
do
    ./mk_image.py $xyz_file
done