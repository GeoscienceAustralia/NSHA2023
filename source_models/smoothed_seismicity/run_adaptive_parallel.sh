#PBS -P w84
#PBS -q express
#PBS -l walltime=03:00:00
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l wd

module load geos/3.8.0
#module load hdf5/1.8.10                                                                                 
module load hdf5/1.10.7
module load openmpi/4.1.4
#module unload python3                                                                                  
#module unload python                                                                                   
module load python3/3.9.2
# Load gdal after python to avoid conflict                                                              
module load gdal/3.5.0

# Local pythonpaths
export PYTHONPATH=.:/home/547/jdg547/.local/lib/python3.9/site-packages:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/modelling/oq-engine/:${PYTHONPATH}

mpirun -np 16 python3 adaptive_smoothing_parallel.py > adaptive_smoothing_parallel.log
