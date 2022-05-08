#!/bin/bash

# ==============================================================================
# Arguments:
# 1 - "sequential", "task-parallel", "data-parallel" or "distributed"
#     depending on job type
# 2, 3, 4 - arguments for qrun
# 5 - path for input file
# 6 - path for results folder
# 7 - number of threads per process (optional)
# 8 - "sort" or "no-sort" to sort edges by weights (optional)
# ==============================================================================

# Adds backslash before each slash
function bs ()
{
  printf "${1//\//\\/}"
}

# ==============================================================================
# Set variables
# ==============================================================================
path_scripts="scripts/"
path_solvers="solvers/"
solver="${path_solvers}${1}.out"
outfile="${6}/$(date +%s)$(ls -1q $6 | wc -l).csv"

if [ $# -lt 7 ]; then
  job_prefix=""
else
  # for "serial" job
  job_prefix="OMP_NUM_THREADS=${7}"
  # for "parallel" job
  export OMP_NUM_THREADS=${7}
fi

if [ $# -lt 8 ];then
  sort=""
else
  sort="$8"
fi

# ==============================================================================
# Prepare job script
# ==============================================================================
# If job type is "serial" (sequential or OpenMP solver)
if [ "$1" = "sequential" ] || [ "$1" = "task-parallel" ] ||
   [ "$1" = "data-parallel"  ]; then

  script="${path_scripts}serial_job.sh"

  # Substitute "$outfile" for stderr output of the SGE script 
  sed "22s/.*/\#\$ -e $(bs $outfile)/" $script > job_tmp.sh

  # Substitute "$infile" for solver input file argument
  sed -i "32s/.*/time ${job_prefix} .\/$(bs $solver) $(bs $5) "$sort"/" job_tmp.sh

# If job type is "parallel" (Open MPI solver)
elif [ "$1" = "distributed" ]; then

  script="${path_scripts}parallel_job.sh"

  # Substitute "$outfile" for stderr output of the SGE script 
  sed "35s/.*/\#\$ -e $(bs $outfile)/" $script > job_tmp.sh

  # Substitute "$infile" for solver input file argument
  sed -i "67s/.*/MY_PARALLEL_PROGRAM=\"$(bs $solver) $(bs $5) "$sort"\"/" job_tmp.sh

else
  exit 1
fi

# ==============================================================================
# Submit job and write result
# ==============================================================================
# Write parameters to outfile
printf "$(date +"%m/%d/%Y %H:%M"), $1, $2, $3, $4, $5, $7, $8, " >> $outfile  

# Submit script
qrun "$2" "$3" "$4" job_tmp.sh

# Clear
rm -f job_tmp.sh

