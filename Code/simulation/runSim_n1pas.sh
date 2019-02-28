### Simple bash shell script to run each all the individual simulation settings

#!/bin/bash

method=seqR
today=`date +%Y-%m-%d.%H:%M:%S`
echo "$today"

## for d in 20 100 500 1000 2000 5000 10000
for d in 20501
do
    R_hat="test_corr_matrix_d="$d"_stand.txt"
    distributions="poisson_data_brca_d="$d"_n=1093.txt"
    cd ~/Research/gpusim/sim_files/slurm_seqR
    mkdir -p $d
    
    for n in 100 1000 2000 5000 10000
    do

	for ((i=1;i<2;i+=1));
	do
	    rep=${i}
	    id=$(echo $n-$d-$method-$rep)
	    echo "working on "$id
	echo "#!/bin/bash

#!/bin/bash	
#SBATCH -n 1
#SBATCH --job-name="$id"
###SBATCH --mem=4GB
#SBATCH --output="$d"/"$id"-"$today".txt

cd ~/Research/gpusim/sim_files

time Rscript ./slurm_seqR/run_seqR.R ./R_hat/"$R_hat" ./distributions/"$distributions" "$d" "$n" "$n$d$rep" ./output_"$method"/"$id"-"$today".txt" > "launch-"$id".slm"
		chmod a+rx "launch-"$id".slm"
		sbatch < "launch-"$id".slm"
	        rm "launch-"$id".slm"
	done
    done
done




