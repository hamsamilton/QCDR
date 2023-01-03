#!/bin/zsh

# Test the current version using bat h 11, using script as a background, with custom inputs and hist and gc inputs
python3 pyRetroPlotter_main.py -ip data/USER_testData_B11.csv -out testingNov12.pdf -gc data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv -hist data/sample11hist_data.csv -bgd data/SCRIPTretro_masterStatistics_allBatches.csv -ctf data/Man_Input_Cutoff_test.xlsx


