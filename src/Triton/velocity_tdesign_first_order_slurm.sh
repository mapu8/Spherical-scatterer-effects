#!/bin/bash -l
#SBATCH --time=0-80:00:00 --mem-per-cpu=1600
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=403
module load matlab
workfolder="/scratch/work/pajunel2/"
outputfile="saved_outputs/velocity_tdesign_first_order_center_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=20
setting="tdesign_first_order_velocity"
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed); exit(0)"
