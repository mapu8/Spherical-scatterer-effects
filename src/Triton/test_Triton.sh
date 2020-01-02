#!/bin/bash -l
#SBATCH --time=0-00:01:00 --mem-per-cpu=500
#SBATCH -p debug
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=0
module load matlab
workfolder="/scratch/work/pajunel2/"
outputfile=saved_outputs/test_$SLURM_ARRAY_TASK_ID.mat
N_diffuse_repetitions=1
setting="2D_vbap"
touch $outputfile
rng_seed=$SLURM_ARRAY_TASK_ID
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$workfolder','$outputfile',$N_diffuse_repetitions,'$setting',$rng_seed); exit(0)"
