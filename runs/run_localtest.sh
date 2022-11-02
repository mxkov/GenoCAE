#!/bin/bash

singularity exec --nv singularity.sif python3 run_gcae.py train \
--datadir data/HOpheno_249ind_1000snp/ --data HOpheno_249ind_1000snp \
--trainedmodeldir ae_out/LOCALTEST/ \
--model_id M1 --train_opts_id ex3 --data_opts_id b_0_4 \
--epochs 10 --save_interval 2 \
2>&1 | tee runs/logs/localtest_train.log

singularity exec --nv singularity.sif python3 run_gcae.py project \
--datadir data/HOpheno_249ind_1000snp/ --data HOpheno_249ind_1000snp \
--trainedmodeldir ae_out/LOCALTEST/ \
--model_id M1 --train_opts_id ex3 --data_opts_id b_0_4 \
--superpops data/HOpheno_249ind_1000snp/HO_superpopulations \
2>&1 | tee runs/logs/localtest_project.log

