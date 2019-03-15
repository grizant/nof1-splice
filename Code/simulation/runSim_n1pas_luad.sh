### Bash shell script to run each all the individual simulation settings
## AG Schissler
## 1 Mar 2019

#!/bin/bash

today=`date +%Y-%m-%d.%H:%M:%S`
echo "$today"
reps=100

## for d in 20 100 500 1000 2000 5000 10000
for dataset in luad
do
    for patientIndex in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58
    do
	for p in 15 30 50 100
	do
	    for pi in 0 0.05 0.1 0.15 0.2
	    do

		id=$(echo $p-$pi-$patientIndex-$dataset-$reps)
		## make unique seed
		datasetID=1
		seed=$datasetID$patientIndex$p$pi
		echo "working on "$id
		echo "#!/bin/bash	
#SBATCH -n 1
#SBATCH --job-name="$id"
###SBATCH --mem=4GB

cd ~/Research/n1pas/sim_files

time Rscript ./runSim_n1pas.R "$dataset" "$patientIndex" "$p" "$pi" "$reps" "$seed" ./output_"$dataset"/"$id"-"$today".csv" > "launch-"$id".slm"
		chmod a+rx "launch-"$id".slm"
		sbatch < "launch-"$id".slm"
	        rm "launch-"$id".slm"
	    done
	done
    done
done
