### Bash shell script to run each all the individual simulation settings
## AG Schissler
## 1 Mar 2019

#!/bin/bash

today=`date +%Y-%m-%d.%H:%M:%S`
echo "$today"
reps=10

## for d in 20 100 500 1000 2000 5000 10000
for dataset in ucec
do
    for patientIndex in 5
    do
	for p in 50
	do
	    for pi in 1
	    do

		id=$(echo $p-$pi-$patientIndex-$dataset-$reps)
		seed=$patientIndex$p$reps$pi
		echo "working on "$id
		echo "#!/bin/bash	
#SBATCH -n 1
#SBATCH --job-name="$id"
###SBATCH --mem=4GB

cd ~/Research/n1pas/sim_files

time Rscript ./runSim_n1pas.R "$dataset" "$patientIndex" "$p" "$pi" "$reps" "$seed" ./output/"$id"-"$today".csv" > "launch-"$id".slm"
		chmod a+rx "launch-"$id".slm"
		## sbatch < "launch-"$id".slm"
	        ##rm "launch-"$id".slm"
	    done
	done
    done
done
