conda activate mdaenv
echo "0"|gmx trjconv -f /home/ricky/Documents/from_work/MD/simulations/production_run//samuli_popc_compare_to_databank/samuli_popc_compare_to_databank.trr -s /home/ricky/Documents/from_work/MD/simulations/production_run//samuli_popc_compare_to_databank/samuli_popc_compare_to_databank.tpr -o samuli_popc_compare_to_databank.xtc -pbc mol -b 00000
./01-calcOrderParameters.py -t  /home/ricky/Documents/from_work/MD/simulations/production_run//samuli_popc_compare_to_databank/samuli_popc_compare_to_databank.gro -x samuli_popc_compare_to_databank.xtc -o samuli_popc_compare_to_databank.dat
rm samuli_popc_compare_to_databank.xtc
awk '$1=="alpha1" {print }' samuli_popc_compare_to_databank.dat.cumulative > samuli_popc_compare_to_databank/alpha1.dat.cumulative
awk '$1=="alpha2" {print }' samuli_popc_compare_to_databank.dat.cumulative > samuli_popc_compare_to_databank/alpha2.dat.cumulative
awk '$1=="beta1" {print }' samuli_popc_compare_to_databank.dat.cumulative > samuli_popc_compare_to_databank/beta1.dat.cumulative
awk '$1=="beta2" {print }' samuli_popc_compare_to_databank.dat.cumulative > samuli_popc_compare_to_databank/beta2.dat.cumulative
rm samuli_popc_compare_to_databank.dat.cumulative
cp samuli_popc_compare_to_databank.dat final_dat_data/samuli_popc_compare_to_databank.dat
rm samuli_popc_compare_to_databank.dat
