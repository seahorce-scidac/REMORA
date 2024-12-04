.. _Perlmutter (NERSC):

Perlmutter (NERSC)
------------------

Recall the GNU Make system is best for use on large computing facility machines and production runs. With the GNU Make implementation, the build system will inspect the machine and use known compiler optimizations explicit to that machine if possible. These explicit settings are kept up-to-date by the AMReX project.

For Perlmutter at NERSC, initialize your environment by sourcing the `saul-env.sh` script in the `Build` directory. For example, from the root of the REMORA directory:

::

   source Build/saul-env.sh

Then follow the general instructions for building REMORA using GNU Make.

   .. note::
      When building, GNU Make may complain that it cannot find NetCDF or MPICH. This is fine.


Building for and running on GPU nodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then build REMORA as, for example (specify your own path to the AMReX submodule in `REMORA/Submodules/AMReX`):

::

   make -j 4 COMP=gnu USE_MPI=TRUE USE_OMP=FALSE USE_CUDA=TRUE AMREX_HOME=$HOME/dev-remora.REMORA/Submodules/AMReX

Finally, you can prepare your SLURM job script, using the following as a guide:

   .. code:: shell

             #!/bin/bash

             ## specify your allocation (with the _g) and that you want GPU nodes
             #SBATCH -A mXXXX_g
             #SBATCH -C gpu

             ## the job will be named "REMORA" in the queue and will save stdout to remora_[job ID].out
             #SBATCH -J REMORA
             #SBATCH -o remora_%j.out

             ## set the max walltime
             #SBATCH -t 10

             ## specify the number of nodes you want
             #SBATCH -N 2

             ## we use the same number of MPI ranks per node as GPUs per node
             #SBATCH --ntasks-per-node=4
             #SBATCH --gpus-per-node=4
             #SBATCH --gpu-bind=none

             # pin to closest NIC to GPU
             export MPICH_OFI_NIC_POLICY=GPU

             # use GPU-aware MPI
             #GPU_AWARE_MPI=""
             GPU_AWARE_MPI="amrex.use_gpu_aware_mpi=1"

             # the -n argument is (--ntasks-per-node) * (-N) = (number of MPI ranks per node) * (number of nodes)
             # set ordering of CUDA visible devices inverse to local task IDs for optimal GPU-aware MPI
             srun -n 8 --cpus-per-task=32 --cpu-bind=cores bash -c "
               export CUDA_VISIBLE_DEVICES=\$((3-SLURM_LOCALID));
               ./REMORA3d.gnu.MPI.CUDA.ex inputs ${GPU_AWARE_MPI}" \
             > test.out

To submit your job script, do `sbatch [your job script]` and you can check its status by doing `squeue -u [your username]`.

Building for and running on CPU nodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then build REMORA as, for example (specify your own path to the AMReX submodule in `REMORA/Submodules/AMReX`):

::

   make -j 4 COMP=gnu USE_MPI=TRUE USE_OMP=TRUE USE_CUDA=FALSE AMREX_HOME=$HOME/dev-remora.REMORA/Submodules/AMReX

Finally, you can prepare your SLURM job script, using the following as a guide:

   .. code:: shell

             #!/bin/bash

             #SBATCH -A mXXXX
             #SBATCH -C cpu
             #SBATCH -q regular

             ## the job will be named "REMORA" in the queue and will save stdout to remora_[job ID].out
             #SBATCH -J REMORA
             #SBATCH -o remora_%j.out

             ## set the max walltime
             #SBATCH -t 10

             ## specify the number of nodes you want
             #SBATCH -N 2

             ## we use 4 ranks per node here as an example. This may not optimize performance
             #SBATCH --ntasks-per-node=4

             ## This configuration assigns one OpenMP thread per physical CPU core.
             ## For this type of thread assignment, we want 128 total threads per node, so we should
             ## have (OMP_NUM_THREADS * ntasks-per-node) = 128
             export OMP_PROC_BIND=spread
             export OMP_PLACES=threads
             export OMP_NUM_THREADS=32

             # the -n argument is (--ntasks-per-node) * (-N) = (number of MPI ranks per node) * (number of nodes)
             srun -n 8 ./REMORA3d.gnu.x86-milan.MPI.OMP.ex inputs > test.out

To submit your job script, do `sbatch [your job script]` and you can check its status by doing `squeue -u [your username]`.
