#!/bin/bash -l
#SBATCH --time=0-120:00:00 --mem-per-cpu=550
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=251,257,272
module load matlab
workfolder="/scratch/work/pajunel2/"
outputfile="saved_outputs/diffuse_tdesign_first_order_center_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=10
setting="tdesign_first_order_diffuse"
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed); exit(0)"
