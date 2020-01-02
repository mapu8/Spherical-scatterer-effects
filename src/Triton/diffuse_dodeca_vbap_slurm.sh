#!/bin/bash -l
#SBATCH --time=0-120:00:00 --mem-per-cpu=550
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=100-199
module load matlab
workfolder="/scratch/work/pajunel2/"
outputfile="saved_outputs/diffuse_dodeca_vbap_d02_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=10
setting="dodeca_vbap_diffuse"
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed); exit(0)"
