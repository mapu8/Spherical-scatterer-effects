#!/bin/bash -l
#SBATCH --time=0-120:00:00 --mem-per-cpu=550
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=400-499
module load matlab
workfolder="/scratch/work/pajunel2/"
outputfile="saved_outputs/diffuse_listening_room_allrad2_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=10
setting="listening_room_allrad_diffuse"
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed); exit(0)"
