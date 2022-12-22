#PBS -P w84
#PBS -q express
#PBS -l walltime=16:00:00
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l wd

#module load intel-cc/12.1.9.293
#module load intel-fc/12.1.9.293
module load geos/3.8.0
#module load hdf5/1.8.10
module load hdf5/1.10.7
module load openmpi/4.1.4
module load gdal
module unload python3
module load python3/3.9.2
#module load python/2.7.11-matplotlib

# To get rtree to run
#export SPATIALINDEX_C_LIBRARY=/short/n74/jdg547/spatialindex-src-1.8.5/lib/libspatialindex_c.so.4
#export LD_LIBRARY_PATH=/short/n74/jdg547/spatialindex-src-1.8.5/lib:$LD_LIBRARY_PATH
# Python paths for local openquake installs and dependencies
export PYTHONPATH=.:/home/547/jdg547/.local/lib/python3.9/site-packages:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/modelling/oq-engine/:${PYTHONPATH}

mpirun -np 16 python3 fixed_smoothing_parallel.py > fixed_smoothing_parallel.log
