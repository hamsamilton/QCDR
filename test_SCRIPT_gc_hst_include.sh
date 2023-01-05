#!/bin/zsh

#Run all script data, with all script data as a background, gc and hist included
python3 pyRetroPlotter_main.py -ip data/SCRIPTretro_masterStatistics_allBatches.csv -out test_script_all.pdf -gc data/SCRIPT_GC_info.csv -hist data/SCRIPT_hist.25_info.csv 

