#!/bin/bash -ue
nohup vmstat 1 10000 > vmstats_SAVER_100_2.txt &
Rscript /home/carisdak/deepimpute_paper_scripts/speed_memory/imputation_runner.R /home/carisdak/deepimpute_paper_scripts/speed_memory/../paper_data/speed_and_memory/mouse1M_100_nonTransposed.csv SAVER 100 2 20
pkill vmstat
