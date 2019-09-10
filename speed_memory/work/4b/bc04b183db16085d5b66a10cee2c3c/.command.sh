#!/bin/bash -ue
nohup vmstat 1 10000 > vmstats_magic_500_2.txt &
python3 /home/carisdak/deepimpute_paper_scripts/speed_memory/imputation_runner.py          --path /home/carisdak/deepimpute_paper_scripts/speed_memory/../paper_data/speed_and_memory/mouse1M_500_transposed.csv          --method magic          --threads 20          --trial 2          --ncells 500
pkill vmstat
