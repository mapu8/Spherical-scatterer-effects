#!/bin/bash -l
#SBATCH --time=0-120:00:00 --mem-per-cpu=550
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=300-399
module load matlab
workfolder="/scratch/work/pajunel2/"
outputfile="saved_outputs/diffuse_tdesign_vbap_d02_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=10
setting="tdesign_vbap_diffuse"
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed); exit(0)"
