#!/bin/bash -ue
cat time* > summary.csv
for f in `ls vmstat*`; do python /home/carisdak/deepimpute_paper_scripts/speed_memory/extract_memory.py $f ; done
