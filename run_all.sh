#!/bin/bash

mkdir -p hrtem_images

for i in {1..99}
do
    process_id=$(printf "%02d" $i)
    mkdir -p tmp_$process_id
    nohup ./run.sh $process_id &
done

cat tmp_*/data.dat > all_data.dat
rm -r tmp_*
