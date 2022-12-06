#!/bin/bash
#SBATCH --gpus 1
#SBATCH -t 1:00:00
#SBATCH --mail-user=max.kovalenko@it.uu.se
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT

# Train with phenomodel

singularity exec --nv singularity.sif python3 run_gcae.py train \
--datadir data/example_tiny/ --data HumanOrigins249_tiny \
--model_id M1 --train_opts_id ex3 --data_opts_id b_0_4 \
--epochs 9999 --patience 300 --save_interval 50 \
--pheno_model_id ph1

singularity exec --nv singularity.sif python3 run_gcae.py project \
--datadir data/example_tiny/ --data HumanOrigins249_tiny \
--model_id M1 --train_opts_id ex3 --data_opts_id b_0_4 \
--superpops data/example_tiny/HO_superpopulations \
--pheno_model_id ph1


# Train without phenomodel

singularity exec --nv singularity.sif python3 run_gcae.py train \
--datadir data/example_tiny/ --data HumanOrigins249_tiny \
--model_id M1 --train_opts_id ex3 --data_opts_id b_0_4 \
--epochs 9999 --patience 300 --save_interval 50

singularity exec --nv singularity.sif python3 run_gcae.py project \
--datadir data/example_tiny/ --data HumanOrigins249_tiny \
--model_id M1 --train_opts_id ex3 --data_opts_id b_0_4 \
--superpops data/example_tiny/HO_superpopulations
