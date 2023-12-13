#!/bin/bash

#$ -m ea
#$ -M dbriers@bu.edu

# This script launches one instance of TRF from an SGE job array
module load anaconda3

# curent directory will be the output directory of TRF!
#cd /restricted/projectnb/johnsonlab/challenge_project_2014/sequence_data/pathoqc_output/split1_output
#cd /scratch
cd /usr4/spclpgm/dbriers/repos/mass/
#SPLIT_DIR='/restricted/projectnb/johnsonlab/challenge_project_2014/sequence_data/pathoqc_output/split1_output'
#split_R1_2.txt

#FASTAFILES=($(cat /restricted/projectnb/johnsonlab/challenge_project_2014/sequence_data/pathoqc_output/split_R1_4.txt))

# Variable $SGE_TASK_ID is a 1-indexed variable unique to each task in the job array
# Bash arrays are 0-indexed, so to get the file we need, must subtract 1
CURR_IDX=${SGE_TASK_ID}

# Use the index calculated above as an index into the FASTAFILES array.
# The array only has the file name, so you need to add the full path.
#CURR_FILE="${FASTAFILES[${CURR_IDX}]}"

# Add the parameters you need before ${CURR_FILE}
# Use quotes around ${CURR_FILE} for safety if filenames have spaces.
# But anyway, try to avoid spaces in file/directory names.

#reset timer
SECONDS=0
#STDERR="/usr4/spclpgm/dbriers/mass_simulations/log_${CURR_IDX}.err"
MASTERLOG="/usr4/spclpgm/dbriers/mass_simulations/master.log"

#RUN TRF and print stdout/err to /dev/null
LABEL='stateMachine8count_7tet_lowerD_numb_early2_002_mid_mc_simple'
SM='state_machines/stateMachine8count_7tet_lowerD_numb_early2_002_mid_mc_simple.csv' #state machine
STDOUT="/usr4/spclpgm/dbriers/mass_simulations/log_${LABEL}_${CURR_IDX}.log"
python /usr4/spclpgm/dbriers/repos/mass/run_simulation.py "${LABEL}_${CURR_IDX}" "${SM}" 1> $STDOUT 2>&1

#mv *.mask ${SPLIT_DIR}

#Track completion time
#F=$(python -c "print $CURR_IDX % 10")
echo "Finished converting ${LABEL}_${CURR_IDX} in $SECONDS seconds." >> ${MASTERLOG}

#if [ $CURR_IDX -eq 1 ]; then
#  jobcomplete "Finished converting ${CURR_IDX} files for job array ${LABEL}."
#fi
