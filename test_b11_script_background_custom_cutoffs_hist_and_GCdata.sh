#!/bin/zsh

# Test the current version using bat h 11, using script as a background, with custom inputs and hist and gc inputs


# run the pipeline and generate the output

# SCRIPT with custom cutoffs
python3 pyRetroPlotter_main.py -ip data/USER_testData_B11.csv -out SCRIPTB11_mancutoffs.pdf -gc data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv -hist data/sample11hist_data.csv -bgd data/SCRIPTretro_masterStatistics_allBatches.csv -ctf data/Man_Input_Cutoff_test.xlsx

# SCRIPT with default cutoffs
python3 pyRetroPlotter_main.py -ip data/USER_testData_B11.csv -out SCRIPTB11_defaultcutoffs.pdf -gc data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv -hist data/sample11hist_data.csv -bgd data/SCRIPTretro_masterStatistics_allBatches.csv 

# SCRIPT with 6 panels
python3 pyRetroPlotter_main.py -ip data/USER_testData_B11.csv -out SCRIPTB11_defaultcutoffs_6panel.pdf -bgd data/SCRIPTretro_masterStatistics_allBatches.csv 

# fibrosis with fibrosis background
python3 pyRetroPlotter_main.py -ip data/fib4pyRetroplotter.csv -out fib_with_fib_background.pdf -gc data/fibGBC_output.csv -hist data/fibhist.254pltr.csv -bgd data/fib4pyRetroplotter.csv

#fibrosis with SCRIPT background
python3 pyRetroPlotter_main.py -ip data/fib4pyRetroplotter.csv -out fib_with_SCRIPT_background.pdf -gc data/fibGBC_output.csv -hist data/fibhist.254pltr.csv -bgd data/SCRIPTretro_masterStatistics_allBatches.csv
# create the folder to store information

FOLDERNAME=all_test_outputs
mkdir $FOLDERNAME

# store the scripts

cp pyRetroPlotter_main.py $FOLDERNAME/pyRetroPlotter_main.py
cp helper_retroFunctions.py $FOLDERNAME/helper_retroFunctions.py

# get version number
cp version_num.txt $FOLDERNAME/version_num.txt

# copy inputs
cp data/USER_testData_B11.csv $FOLDERNAME/USER_testData_B11.csv
cp data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv $FOLDERNAME/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv
cp data/sample11hist_data.csv $FOLDERNAME/sample11hist_data.csv
cp data/SCRIPTretro_masterStatistics_allBatches.csv $FOLDERNAME/SCRIPTretro_masterStatistics_allBatches.csv
cp data/Man_Input_Cutoff_test.xlsx $FOLDERNAME/Man_Input_Cutoff_test.xlsx
cp data/fib4pyRetroplotter.csv $FOLDERNAME/fib4pyRetroplotter.csv
cp data/fibGBC_output.csv $FOLDERNAME/fibGBC_output.csv
cp data/fibhist.254pltr.csv $FOLDERNAME/fibhist_25_4pltr.csv

# copy outputs

cp SCRIPTB11_mancutoffs.pdf $FOLDERNAME/SCRIPTB11_mancutoffs.pdf
cp SCRIPTB11_defaultcutoffs.pdf $FOLDERNAME/SCRIPTB11_defaultcutoffs.pdf
cp SCRIPTB11_defaultcutoffs_6panel.pdf $FOLDERNAME/SCRIPTB11_defaultcutoffs_6panel.pdf
cp fib_with_fib_background.pdf $FOLDERNAME/fib_with_fib_background.pdf
cp fib_with_SCRIPT_background.pdf $FOLDERNAME/fib_with_SCRIPT_background.pdf






