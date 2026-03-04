#!/bin/bash

export path_xyz="/home/victor.barrere@crmd.cnrs-orleans.fr/Documents/Data/Data_xyz/"
export path_processed="/home/victor.barrere@crmd.cnrs-orleans.fr/Documents/Data/Data_processed/"
export path_new_xyz=$path_processed"XYZ/"
path_hrtem='/home/victor.barrere@crmd.cnrs-orleans.fr/Documents/Data/Data_HRTEM'
path_msa_prm=$path_hrtem/msa_prm
path_wavimg_prm=$path_hrtem/wavimg_prm
path_img=$path_hrtem/HRTEM_image

# (Devient new path) path_xyz='/home/victor.barrere@crmd.cnrs-orleans.fr/Documents/Data/Data_processed/XYZ'

# Creation of the new XYZ files with the cluster separation
mkdir -p $path_new_xyz
rm -f $path_processed/data.dat

$ovitos cluster_separation.py
$ovitos particle_analyze.py


# Image creation
export electron_energy=300
export nx=64
export ny=64
export nz=10
export Nlp=0
export Nabr=0

mkdir -p $path_hrtem

for file_xyz in $path_new_xyz/*.xyz
do

    export file_xyz=$file_xyz
    export file_cel=${file_xyz%.xyz}.cel

    id_sim=$(basename "${file_xyz%.xyz}")
    path_id_sim=$path_hrtem/$id_sim

    export file_sli=$path_id_sim/slice
    export file_msa_prm=$path_msa_prm/msa_$id_sim.prm
    export file_wav=$path_id_sim/HRTEM_image_$id_sim.wav
    export file_wav_prm=$path_wavimg_prm/wavimg_$id_sim.prm
    export file_image=$path_id_sim/HRTEM_image_$id_sim.dat
    export final_img=$path_img/$id_sim.png

    mkdir -p $path_id_sim $path_msa_prm $path_wavimg_prm $path_img

    ./xyz_to_cel.py
    ./create_msa_prm.py
    ./create_wav_prm.py
    
    celslc -cel $file_cel -slc $file_sli -nx $nx -ny $ny -nz $nz -ht $electron_energy -dwf -abs
    rm $path_id_sim/slice.prm
    
    msa -prm $file_msa_prm -out $file_wav /ctem
    rm $path_id_sim/*.sli $path_new_xyz/*.cel
    mv  $path_id_sim/HRTEM_image_*.wav $file_wav
    
    wavimg -prm $file_wav_prm -out $file_image
    
    ./hrtem_image.py
    rm -r $path_id_sim

done