#!/bin/bash
#SBATCH -A CSC143
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --core-spec=0
#SBATCH -J CPUJHTDB
#SBATCH -o CPUJHTDB.out

# Suppose you've successfully downloaded and sliced data into JHTDB (1024x2048x2048, [1536:2560, 1024:3072, 1024:3072] from 4096x4096x4096 isotropic4096 data)
# You have 64 CPUs
# Slice the JHTDB VelocityX, VelocityY, VelocityZ into 64 pieces 256x512x512
# You have a directory named JHTDB with in it exists under current directory.
# make sure you have enough space (~128GB) to store JHTDB and refactored data

set -x
set -e

a1=0.1
r=0.1
n=5
error_bounds=()

a=$a1
for ((i = 1; i <= n; i++)); do
    error_bounds+=($a)
    a=$(echo "scale=10; $a * $r" | bc)
done

error_bounds=($(printf "%s\n" "${error_bounds[@]}" | sort -nr))

refactor="./build/parallel_src/para_refactor"
reconstructor="./build/parallel_src/para_VTOT"

output_file="CPU_JHTDB_output.txt"
tmp_file="CPU_JHTDB_tmp_output.txt"
>$output_file
>$tmp_file

SRUN="srun -A CSC143 -N 1 -n 64 --ntasks-per-node=64"

$SRUN $refactor 2 JHTDB JHTDB/JHTDB > $tmp_file

time=$(grep "max_elapsed_time" $tmp_file | head -n 1)
echo "Refactor: $time" >> $output_file

for error_bound in "${error_bounds[@]}"; do
    $SRUN $reconstructor $error_bound JHTDB/JHTDB retrieved_data > $tmp_file
    bitrate=$(grep "bitrate" $tmp_file | head -n 1)
    readtime=$(grep "IO_time" $tmp_file | head -n 1)
    time=$(grep "elapsed_time" $tmp_file | head -n 1)
    requested_max_error=$(grep "Target" $tmp_file | head -n 1)
    est_max_error=$(grep "est_error" $tmp_file | head -n 1)
    real_max_error=$(grep "act_error" $tmp_file | head -n 1)
    echo "Request eb = $error_bound, $bitrate, $readtime, $time, $requested_max_error, $est_max_error, $real_max_error" >> $output_file
done

cat $output_file
rm $tmp_file
rm $output_file
