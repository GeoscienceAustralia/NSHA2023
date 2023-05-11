#PBS -P w84
#PBS -q express
#PBS -l walltime=24:00:00
#PBS -l ncpus=1
#PBS -l mem=32GB
#PBS -l wd

module load geos/3.8.0                                                                                                   
module load hdf5/1.10.7
module load openmpi/4.1.4                                                                                                  
module load python3/3.9.2
# Load gdal after python to avoid conflict                                                                                
module load gdal/3.5.0


# Local pythonpaths                                                                                                        
export PYTHONPATH=.:/home/547/jdg547/.local/lib/python3.9/site-packages:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/:${PYTHONPATH}
export PYTHONPATH=.::/home/547/jdg547/modelling/oq-engine/:${PYTHONPATH}

python3 combine_ss_frankel_models.py >& combine_ss_frankel_models.log
#python3 combine_ss_models.py >& combine_ss_models.log         
