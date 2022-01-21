# multiscale-power-system-expansion

This is a repo for the paper "A Nested Cross Decomposition Algorithm for Power System Capacity Expansion with Multiscale Uncertainties".

### Introduction
Modern electric power systems have witnessed rapidly increasing penetration of renewable energy, storage, electrical vehicles and various demand response resources. The electric infrastructure planning is thus facing more challenges due to the variability and uncertainties arising from the diverse new resources. We develop a multistage and multiscale stochastic mixed integer programming model to capture both the coarse-temporal-scale uncertainties, such as investment cost and long-run demand stochasticity, and fine-temporal-scale uncertainties, such as hourly renewable energy output and electricity demand uncertainties, for the power system capacity expansion problem. To be applied to a real power system, the resulting model will lead to extremely large-scale mixed integer programming problems, which suffer not only the well-known curse of dimensionality, but also computational difficulties with a vast number of integer variables at each stage. In addressing such challenges associated with the developed model, we propose a nested cross decomposition (NCD) algorithm that consists of two layers of decomposition, that is, the Dantzig-Wolfe decomposition and L-shaped decomposition. 

### The repo includes:
* the C++ MPI codes programmed for the NCD algorithm
* the data for a tested 6-bus power system 
* the data for a modified IEEE 118-bus pow system
* the instances used in the numerical experiments

### Folder directory
#### Solver
* the C++ MPI codes programmed for the NCD algorithm
#### Instances
* the data for a tested 6-bus power system and corresponding numerical results 
* the data for a modified IEEE 118-bus pow system and corresponding numerical results 
#### Results
* More numerical results for all the test instances on the IEEE 118-bus system

### To compile
One can use CmakeLists.txt, Makefile or the solution file contained in the "Solver" folder to compile the program. 

### To Run
* serial computing: ./Model directory_to_test_instance_input_folder -n
* parallel computing: mpiexec -n number_of_ranks Model directory_to_test_instance_input_folder -n
* a slurm script should be used if the parallel computing job is submitted to a HPC cluster with the slurm workload manager