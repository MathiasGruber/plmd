## Note: Created from runAll.py

#!/bin/sh
#
#PBS -N cases/5
#PBS -l nodes=1:ppn=8
#PBS -l walltime=72:00:00

cd $PBS_O_WORKDIR

# Set absolute path to amber
amber=$AMBERHOME/bin/sander.MPI

# Load mpi
module load mpi/gcc-4.7.2-openmpi-1.6.3

# Number of MD steps
md_steps=50

######################

# Minimize the system
echo "Minimization 1; ligand & peptide restrained"
if [[ -e cases/5/md-files/min.rst ]]; then
    echo '-- min.rst file exists, skipping minimization 1'
else
    echo '-- Running Minimization 1'
    mpirun $amber -O -i cases/5/in_files/explicit_min.in -o cases/5/md-logs/outMin.log -c cases/5/md-files/peptide.inpcrd -p cases/5/md-files/peptide.prmtop -r cases/5/md-files/min.rst -ref cases/5/md-files/peptide.inpcrd
    
    # Create pdb file for the minimized structure
    ${AMBERHOME}/bin/ambpdb -p cases/5/md-files/peptide.prmtop < cases/5/md-files/min.rst > cases/5/pdb-files/minimization.pdb
fi

######################

echo "Heating, no restraints"
if [[ -e cases/5/md-files/equil0.rst ]]; then
    echo '-- equil0.rst file exists, skipping heating'
elif [[ -e cases/5/md-files/min.rst ]]; then
    echo '-- Running Heating'
    mpirun $amber -O -i cases/5/in_files/explicit_heat.in -o cases/5/md-logs/outHeat.log -c cases/5/md-files/min.rst -p cases/5/md-files/peptide.prmtop -r cases/5/md-files/equil0.rst -x cases/5/md-files/heat.mdcrd -ref cases/5/md-files/min.rst
    
    # Convert the resulting structure to a pdb-file
    ${AMBERHOME}/bin/ambpdb -p cases/5/md-files/peptide.prmtop < cases/5/md-files/equil0.rst > cases/5/pdb-files/heating.pdb
else
    echo '-- min.rst not found, skipping heating'
fi

######################

# A function for running another step of the script
runMD(){

# define the local step nr, previous step nr, and new step nr
local cStep=$1
local pStep=$(($1-1))

# Start amber MD simulation
mpirun $amber -O -i cases/5/in_files/explicit_equil.in -o cases/5/md-logs/outMD${cStep}.log -c cases/5/md-files/equil${pStep}.rst -p cases/5/md-files/peptide.prmtop -r cases/5/md-files/equil${cStep}.rst -x cases/5/md-files/equil${cStep}.mdcrd

# Create pdb file from results
$AMBERHOME/bin/ambpdb -p cases/5/md-files/peptide.prmtop < cases/5/md-files/equil${cStep}.rst > cases/5/pdb-files/equil${cStep}.pdb
}

# Run the MD simulations. Figure out last stop, and take "md_steps" runs from that
arr=( md-files/equil*.rst )  # * is list of all file and dir names
n=${#arr[@]}
echo "Molecular Dynamics"
echo "-- Number of previous MD runs: "$(($n-1))

startI=$n
endI=$(($n+$md_steps))

while [ $startI -lt $endI ]
do
    if [[ -e cases/5/md-files/equil$(($startI+1)).rst ]]; then
        echo "-- #"$startI" MD already run, skipping to next"
    elif [[ -e cases/5/md-files/equil$(($startI-1)).rst ]]; then
        echo "-- Running script #"$startI
        runMD $startI
    else
        echo '-- previous equil file not found, skipping equilibration #'$startI
    fi
    startI=`expr $startI + 1`
done