#!/bin/bash -l
#SBATCH --time=0-80:00:00 --mem-per-cpu=9000
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=700
module load matlab
workfolder="/scratch/work/pajunel2/"
outputfile="saved_outputs/vbap_2D_16speakers_velocity_noscat_full_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=20
setting="2D_vbap_16_velocity_noscat"
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed); exit(0)"
