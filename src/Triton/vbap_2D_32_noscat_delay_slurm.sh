#!/bin/bash -l
#SBATCH --time=0-80:00:00 --mem-per-cpu=550
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=523
module load matlab
workfolder="/scratch/work/pajunel2/"
delay_distance=0.02
outputfile="saved_outputs/vbap_2D_32speakers2_noscat_delay_"$delay_distance"_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=20
setting="2D_vbap_32_noscat_delay"
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed,$delay_distance); exit(0)"
