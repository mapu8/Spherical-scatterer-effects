#!/bin/bash -l
#SBATCH --time=0-00:05:00 --mem-per-cpu=800
#SBATCH -p debug
#SBATCH -o ./logs/job-%a.out
#SBATCH --array=0-9
module load matlab
outputfile="saved_outputs/diffuse_reference_$SLURM_ARRAY_TASK_ID.mat"
N_diffuse_repetitions=3
touch $outputfile
matlab -nojvm -r "savePressureFieldsTriton('$outputfile',$N_diffuse_repetitions); exit(0)"
