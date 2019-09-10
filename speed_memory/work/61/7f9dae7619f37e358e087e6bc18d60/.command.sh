#!/bin/bash -ue
nohup vmstat 1 10000 > vmstats_magic_100_1.txt &
python3 /home/carisdak/deepimpute_paper_scripts/speed_memory/imputation_runner.py          --path /home/carisdak/deepimpute_paper_scripts/speed_memory/../paper_data/speed_and_memory/mouse1M_100_transposed.csv          --method magic          --threads 20          --trial 1          --ncells 100
pkill vmstat
