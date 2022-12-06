#!/bin/bash
#SBATCH --gpus 1
#SBATCH -t 1:00:00
#SBATCH --mail-user=max.kovalenko@it.uu.se
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT

datadir=data/example_tiny/
prefix=HumanOrigins249_tiny
superpops=data/example_tiny/HO_superpopulations
mid=M1
tid=ex3
did=b_0_4

phmodels=('ph0' 'ph1')
for phid in ${phmodels[@]}; do
	for k in {1..5}; do
		singularity exec --nv singularity.sif python3 run_gcae.py train \
		--datadir $datadir --data $prefix \
		--model_id $mid --train_opts_id $tid --data_opts_id $did \
		--epochs 9999 --patience 300 --save_interval 50 \
		--pheno_model_id $phid \
		2>&1 | tee $phid.$k.train.log

		singularity exec --nv singularity.sif python3 run_gcae.py project \
		--datadir $datadir --data $prefix \
		--model_id $mid --train_opts_id $tid --data_opts_id $did \
		--superpops $superpops \
		--pheno_model_id $phid \
		2>&1 | tee $phid.$k.project.log

		outdir=ae_out/ae.$mid.$phid.$tid.$did.$prefix
		mv $phid.$k.* $outdir
		mv $outdir $outdir.$k
	done
done

