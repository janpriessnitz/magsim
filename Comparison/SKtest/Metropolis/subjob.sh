#PBS -N Metropolis
#PBS -A OPEN-23-38
#PBS -q qprod
#PBS -l select=1:ncpus=128
#PBS -l walltime=48:00:00
#PBS -j oe

# additional info files
JOB_INFO=$PBS_O_WORKDIR/info.out
STD_OUT=$PBS_O_WORKDIR/uppasd.out

# setting the automatical cleaning of the SCRATCH
# trap 'clean_scratch' TERM EXIT

# load modules
module load iccifort/2020.4.304
module load impi/2019.9.304-iccifort-2020.4.304
module load imkl/2021.2.0-iimpi-2021a

# create scratch
SCR="/lscratch/$PBS_JOBID"
mkdir -p $SCR

# copy input files to the scratch
cp $PBS_O_WORKDIR/inpsd.dat $SCR
cp $PBS_O_WORKDIR/jfile $SCR
cp $PBS_O_WORKDIR/dmfile $SCR
cp $PBS_O_WORKDIR/kfile $SCR
cp $PBS_O_WORKDIR/posfile $SCR
cp $PBS_O_WORKDIR/momfile $SCR

# change directory to the sratch
cd $SCR || exit

# run simulation
/home/balaz/UppASD/UppASD-master/source/sd

# archive data
archive="results.tar"
tar -cvf $archive *.out
gzip $archive

# copy results to my Datadir
cp ./$archive".gz" $PBS_O_WORKDIR

# remove files from the scratch
rm *

