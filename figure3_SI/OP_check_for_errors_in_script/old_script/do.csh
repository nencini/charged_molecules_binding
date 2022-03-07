#rermove old scripts


rm calculate_OP.csh
mkdir final_dat_data

directory='/home/ricky/Documents/from_work/MD/simulations/production_run/'

#make the script to analyze the OP change with ion concentration
while IFS="," read -r f1
do
	mkdir $f1
	
       echo 'conda activate mdaenv' >> calculate_OP.csh
       echo 'echo "0"|gmx trjconv -f '$directory'/'$f1'/'$f1'.xtc -s '$directory'/'$f1'/'$f1'.tpr -o '$f1'.xtc -pbc mol -b 200000' >> calculate_OP.csh 
       echo './01-calcOrderParameters.py -t  '$directory'/'$f1'/'$f1'.gro -x '$f1'.xtc -o '$f1'.dat' >> calculate_OP.csh
       echo 'rm '$f1'.xtc' >> calculate_OP.csh

	echo "awk '\$1=="'"'"alpha1"'"'" {print }' "$f1".dat.cumulative > "$f1"/alpha1.dat.cumulative" >> calculate_OP.csh
       echo "awk '\$1=="'"'"alpha2"'"'" {print }' "$f1".dat.cumulative > "$f1"/alpha2.dat.cumulative" >> calculate_OP.csh
       echo "awk '\$1=="'"'"beta1"'"'" {print }' "$f1".dat.cumulative > "$f1"/beta1.dat.cumulative" >> calculate_OP.csh
       echo "awk '\$1=="'"'"beta2"'"'" {print }' "$f1".dat.cumulative > "$f1"/beta2.dat.cumulative" >> calculate_OP.csh

      



       echo 'rm '$f1'.dat.cumulative' >> calculate_OP.csh
       echo 'cp '$f1'.dat final_dat_data/'$f1'.dat' >> calculate_OP.csh
       echo 'rm '$f1'.dat' >> calculate_OP.csh




done < simulations.dat






chmod +x calculate_OP.csh


