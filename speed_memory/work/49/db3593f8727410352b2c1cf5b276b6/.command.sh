#!/bin/bash -ue
nohup vmstat 1 10000 > vmstats_SAVER_500_1.txt &
Rscript /home/carisdak/deepimpute_paper_scripts/speed_memory/imputation_runner.R /home/carisdak/deepimpute_paper_scripts/speed_memory/../paper_data/speed_and_memory/mouse1M_500_nonTransposed.csv SAVER 500 1 20
pkill vmstat
