#!/bin/bash
#SBATCH -A C3SE2023-1-17 -p vera
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=gustajo@chalmers.se
#SBATCH --mail-type=end  # send mail when job ends


# instead of using the "SBATCH -o run_init-X.log" line here,
# include the flag "-o run_init-%A_%a.log" in the sbatch submission command:
# sbatch -o logs/run_evaless_depmap-%A-%a.log --array=1-2 run_evaless_depmap.sh
#change to 1-40 later
module load MATLAB/2019a
module load Gurobi/9.5.2

#use the upper line the first time it is run. If it generates some models and then for example runs out of time, the second line can be used to continue on the previous run
matlab -nodesktop -nodisplay -r "getTaskEssentialGenesCluster('./data/depmap_models-', ${SLURM_ARRAY_TASK_ID}); exit" < /dev/null &

wait

