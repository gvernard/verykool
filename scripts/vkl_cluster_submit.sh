#!/bin/bash
 
#SBATCH --partition=compute                             # Partition to use 
#SBATCH --nodes=2                                       # Number of nodes to request 
#SBATCH --tasks-per-node=48                             # Number of cores per node
#SBATCH --time=10:00:00                                 # Format is DD-HH:MM:SS
#SBATCH --mem=20480MB                                   # Memory per node. Default is 4 GB for on-prem nodes
                                                        # and aws-t2 queues. Default is 2 GB for aws-c5.4x and
                                                        # aws.c5-12x queues. Values are specified in MBs.
                                                        # Possible to also use Can use K,G,T. Setting to 0 will
                                                        # request all the memory on a node
#SBATCH --output=%x-%J.out                              # Name of file to send standard output to. You should
                                                        # also send output from your programs using their
                                                        # output options if available
#SBATCH --error=%x-%J.err                               # Name of file to send standard error to. You should
                                                        # also send errors from your programs using their
                                                        # error output options if available
PBS_NODEFILE=$(/usr/bin/generate_pbs_nodefile)          # Generate hostfile for use with mpirun command
NPROCS=$(wc -l < $PBS_NODEFILE)                         # Get total number of cores for mpirun command
MYPATH=$1
RUN=$2

module load OpenMPI/openmpi-4.0.4
module load GNU/gnu-gsl-2.7.1
module load Boost/boost_1_75_0

yes | ${HOME}/verykool/bin/vkl_driver ${MYPATH} ${RUN} ${PBS_NODEFILE} ${NPROCS}
#yes | ${HOME}/verykool/bin/vkl_driver /home/gvernardos/old_s22_unfiltered/ curv_evi/ ${PBS_NODEFILE} ${NPROCS}
