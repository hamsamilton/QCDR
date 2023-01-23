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

# save the custom cutoffs script
CUSTNAME=$FOLDERNAME/SCRIPTB11_mancutoffs
INPUTFOLDER=$CUSTNAME/inputs
OUTPUTFOLDER=$CUSTNAME/outputs
echo $CUSTNAME
echo $INPUTFOLDER
echo $OUTPUTFOLDER

mkdir $CUSTNAME
mkdir $INPUTFOLDER
mkdir $OUTPUTFOLDER
cp data/USER_testData_B11.csv $INPUTFOLDER/USER_testData_B11.csv
cp data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv $INPUTFOLDER/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv
cp data/sample11hist_data.csv $INPUTFOLDER/sample11hist_data.csv
cp data/SCRIPTretro_masterStatistics_allBatches.csv $INPUTFOLDER/SCRIPTretro_masterStatistics_allBatches.csv
cp data/Man_Input_Cutoff_test.xlsx $INPUTFOLDER/Man_Input_Cutoff_test.xlsx
cp SCRIPTB11_mancutoffs.pdf $OUTPUTFOLDER/SCRIPTB11_mancutoffs.pdf
echo "python3 pyRetroPlotter_main.py -ip data/USER_testData_B11.csv -out SCRIPTB11_mancutoffs.pdf -gc data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv -hist data/sample11hist_data.csv -bgd data/SCRIPTretro_masterStatistics_allBatches.csv -ctf data/Man_Input_Cutoff_test.xlsx" > $CUSTNAME/pythoncall.txt 

# save default cutoffs SCRIPT
CUSTNAME=$FOLDERNAME/SCRIPTB11_defaultcutoffs
INPUTFOLDER=$CUSTNAME/inputs
OUTPUTFOLDER=$CUSTNAME/outputs
mkdir $CUSTNAME
mkdir $INPUTFOLDER
mkdir $OUTPUTFOLDER
cp data/USER_testData_B11.csv $INPUTFOLDER/USER_testData_B11.csv
cp data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv $INPUTFOLDER/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv
cp data/sample11hist_data.csv $INPUTFOLDER/sample11hist_data.csv
cp data/SCRIPTretro_masterStatistics_allBatches.csv $INPUTFOLDER/SCRIPTretro_masterStatistics_allBatches.csv
cp SCRIPTB11_defaultcutoffs.pdf $OUTPUTFOLDER/SCRIPTB11_defaultcutoffs.pdf
echo "python3 pyRetroPlotter_main.py -ip data/USER_testData_B11.csv -out SCRIPTB11_defaultcutoffs.pdf -gc data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv -hist data/sample11hist_data.csv -bgd data/SCRIPTretro_masterStatistics_allBatches.csv" > $CUSTNAME/python_call.txt
#


# save 6panel test SCRIPT
CUSTNAME=$FOLDERNAME/SCRIPTB11_defaultcutoffs_6panel
INPUTFOLDER=$CUSTNAME/inputs
OUTPUTFOLDER=$CUSTNAME/outputs
mkdir $CUSTNAME
mkdir $INPUTFOLDER
mkdir $OUTPUTFOLDER
cp data/USER_testData_B11.csv $INPUTFOLDER/USER_testData_B11.csv
cp data/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv $INPUTFOLDER/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv
cp data/sample11hist_data.csv $INPUTFOLDER/sample11hist_data.csv
cp data/SCRIPTretro_masterStatistics_allBatches.csv $INPUTFOLDER/SCRIPTretro_masterStatistics_allBatches.csv
cp SCRIPTB11_defaultcutoffs_6panel.pdf $OUTPUTFOLDER/SCRIPTB11_defaultcutoffs_6panel.pdf

echo "python3 pyRetroPlotter_main.py -ip data/USER_testData_B11.csv -out SCRIPTB11_defaultcutoffs_6panel.pdf -bgd data/SCRIPTretro_masterStatistics_allBatches.csv" > $CUSTNAME/pythoncall.txt 

# save fibrosis default background
CUSTNAME=$FOLDERNAME/fib_with_fib_background
INPUTFOLDER=$CUSTNAME/inputs
OUTPUTFOLDER=$CUSTNAME/outputs
mkdir $CUSTNAME
mkdir $INPUTFOLDER
mkdir $OUTPUTFOLDER
cp data/fib4pyRetroplotter.csv $INPUTFOLDER/fib4pyRetroplotter.csv
cp data/fibGBC_output.csv $INPUTFOLDER/fibGBC_output.csv
cp data/fibhist.254pltr.csv $INPUTFOLDER/fibhist_25_4pltr.csv
cp fib_with_fib_background.pdf $OUTPUTFOLDER/fib_with_fib_background.pdf
echo "python3 pyRetroPlotter_main.py -ip data/fib4pyRetroplotter.csv -out fib_with_fib_background.pdf -gc data/fibGBC_output.csv -hist data/fibhist.254pltr.csv -bgd data/fib4pyRetroplotter.csv" > $CUSTNAME/pythoncall.txt 

#save fibrosis with script background
CUSTNAME=$FOLDERNAME/fib_with_SCRIPT_background
INPUTFOLDER=$CUSTNAME/inputs
OUTPUTFOLDER=$CUSTNAME/outputs
mkdir $CUSTNAME
mkdir $INPUTFOLDER
mkdir $OUTPUTFOLDER
cp data/fib4pyRetroplotter.csv $INPUTFOLDER/fib4pyRetroplotter.csv
cp data/SCRIPTretro_masterStatistics_allBatches.csv $INPUTFOLDER/SCRIPTretro_masterStatistics_allBatches.csv
cp data/fibGBC_output.csv $INPUTFOLDER/fibGBC_output.csv
cp data/fibhist.254pltr.csv $INPUTFOLDER/fibhist_25_4pltr.csv
cp fib_with_SCRIPT_background.pdf $OUTPUTFOLDER/fib_with_SCRIPT_background.pdf
echo "python3 pyRetroPlotter_main.py -ip data/fib4pyRetroplotter.csv -out fib_with_SCRIPT_background.pdf -gc data/fibGBC_output.csv -hist data/fibhist.254pltr.csv -bgd data/SCRIPTretro_masterStatistics_allBatches.csv" > $CUSTNAME/pythoncall.txt

SCRIPTFOLDER=$FOLDERNAME/scripts
mkdir $SCRIPTFOLDER
cp pyRetroPlotter_main.py $SCRIPTFOLDER/pyRetroPlotter_main.py
cp helper_retroFunctions.py $SCRIPTFOLDER/helper_retroFunctions.py

# get version number
cp version_num.txt $FOLDERNAME/version_num.txt

# send git log to the folder
git log > $FOLDERNAME/git_log.txt




