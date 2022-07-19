import matplotlib
import matplotlib.gridspec as gridspec
import os
import argparse

matplotlib.use('Agg')
import seaborn
import matplotlib.pyplot as plt
import sqlite3
import sys
import json
import seaborn as sns
import time
import shutil
import glob
import csv
import sqlalchemy as sq
import subprocess
import numpy as np
import pandas as pd
import xlsxwriter
from sklearn.linear_model import LinearRegression
from itertools import zip_longest
from scipy.stats import norm
from scipy import stats
from sklearn.preprocessing import LabelEncoder
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import MinMaxScaler
import helper_retroFunctions
import matplotlib.offsetbox
from datetime import datetime
# from pandas.plotting import register_matplotlib_converters
from statsmodels.stats.weightstats import ztest






_masterDB = "/projects/b1079/tools/U19_masterDB_dev/U19_masterDB.sqlite"
_pipeline_dir = "/projects/b1079/U19_RNAseq_Pipeline_SLURM_v01"

# Define FAIL/WARN thresholds for plots
##FAIL
_ipReads_cutoff_fail = 5
_trimmedReads_cutoff_fail = 80
_uniqAligned_cutoff_fail = 50
_exonMapping_cutoff_fail = 40
_riboScatter_cutoff_fail = 0.50
_violin_cutoff_fail = 1
_violin_cutoff_adapter_fail = 1
_violin_cutoff_overrep_fail = 1

##WARN
_ipReads_cutoff_warn = 10
_trimmedReads_cutoff_warn = 90
_uniqAligned_cutoff_warn = 60
_exonMapping_cutoff_warn = 50
_riboScatter_cutoff_warn = 0.35
_violin_cutoff_warn = 0.5
_violin_cutoff_adapter_warn = 0.5
_violin_cutoff_overrep_warn = 0.5

required_modules = ['python/anaconda']
_notification_email = "gaurav.gadhvi@northwestern.edu"


def join_levels(_df):
    for i, col in enumerate(_df.columns.levels):
        columns = np.where(col.str.contains('Unnamed'), '', col)
        #print(columns)

        _df.columns.set_levels(columns, level=i, inplace=True)

    # print([' '.join(col).strip() for col in _df.columns.values])
    # print(_df.columns)

    _df.columns = [' '.join(_col).strip() for _col in _df.columns.values]

    return _df

def create_connection(_db_file):
    try:
        _conn = sqlite3.connect(_db_file)
        return _conn
    except Error as e:
        print(e)

    return None



def rest_Dev():

    
    '''{BEGIN Current Batch Input from USER'''

    ######## Take this input from user (using specified format/columns) #######

    os.chdir(_projectDir + "/Output/3-STAROutput_Bam/stats")

    _currBatch_stat = pd.DataFrame()

    for _log in glob.glob("*sequencingStats.json"):
        _currBatch_stat = pd.read_json(_log, orient='split')

        _runDate = _log.split("__")[1].split("_hg19")[0]
        _runDate = pd.to_datetime(_runDate).strftime("%m-%d-%Y")

        _currBatch_stat["Date"] = _runDate
        _currBatch_stat["Project"] = str(_projectName)

    _currBatch_stat = _currBatch_stat.loc[:,
                      ['SampleID', '# Sequenced Reads', '% of Reads after Trimming', '# Uniquely Aligned Reads',
                       '% Uniquely Aligned Reads', '% Mapped to Exons/Aligned',
                       '# Uniquely Aligned Reads overlaping rRNA',
                       '# Expressed Genes Detected',
                       'Project', 'Date']]

    _stat_cols = ['Sample', 'Input_Size', 'Percent_PostTrim', 'Num_Uniquely_Aligned', 'Percent_Uniquely_Aligned',
                  'Percent_Exonic', 'Num_Uniquely_Aligned_rRNA', 'Num_expressed_genes', 'Project', 'Date']

    _currBatch_stat.columns = _stat_cols

    _currBatch_stat.loc[:, 'Input_Size'] = pd.to_numeric(_currBatch_stat['Input_Size'])
    _currBatch_stat.loc[:, 'Percent_PostTrim'] = pd.to_numeric(_currBatch_stat['Percent_PostTrim'].str.strip("%"))
    _currBatch_stat.loc[:, 'Num_Uniquely_Aligned'] = pd.to_numeric(_currBatch_stat['Num_Uniquely_Aligned'])
    _currBatch_stat.loc[:, 'Percent_Uniquely_Aligned'] = pd.to_numeric(
        _currBatch_stat['Percent_Uniquely_Aligned'].str.strip("%"))
    _currBatch_stat.loc[:, 'Percent_Exonic'] = pd.to_numeric(_currBatch_stat['Percent_Exonic'].str.strip("%"))
    _currBatch_stat.loc[:, 'Num_Uniquely_Aligned_rRNA'] = pd.to_numeric(_currBatch_stat['Num_Uniquely_Aligned_rRNA'])
    _currBatch_stat.loc[:, 'Num_expressed_genes'] = pd.to_numeric(_currBatch_stat['Num_expressed_genes'])
    _currBatch_stat.loc[:, 'Date'] = pd.to_datetime(_currBatch_stat['Date'], format='%m-%d-%Y').dt.date.astype('O')

    # print(_currBatch_stat)

    ## Read UntrimmedQC stats as a separate DataFrame and merge with _currBatch_stat
    _currBatch_qc_untrim = pd.read_excel(
        _projectDir + "/Output/1-UntrimmedQC/QC_Table/" + _projectName + "_Untrimmed_FastQC_Table.xlsx",
        index_col=0, header=[1, 2])

    # _currBatch_qc_untrim = join_levels(_currBatch_qc_untrim)
    # _currBatch_qc_untrim = join_levels_mat3(_currBatch_qc_untrim)
    _currBatch_qc_untrim = join_levels(_currBatch_qc_untrim)

    # print(_currBatch_qc_untrim.columns)

    _currBatch_qc_untrim = _currBatch_qc_untrim.loc[:,
                           ['Sample', 'Overrepresented Sequences [% of Total Sequences]', '% of Adapter Content']]
    _currBatch_qc_untrim.columns = ['Sample', 'Percent_Overrepresented_Seq_Untrimmed',
                                    'Percent_Adapter_Content_Untrimmed']

    _currBatch_qc_untrim.loc[:, 'Percent_Overrepresented_Seq_Untrimmed'] = pd.to_numeric(
        _currBatch_qc_untrim['Percent_Overrepresented_Seq_Untrimmed'].str.strip("%"))
    _currBatch_qc_untrim.loc[:, 'Percent_Adapter_Content_Untrimmed'] = pd.to_numeric(
        _currBatch_qc_untrim['Percent_Adapter_Content_Untrimmed'].str.strip("%"))

    # print(_currBatch_qc_untrim)

    ## Read TrimmedQC stats as a separate DataFrame and merge with _currBatch_stat
    _currBatch_qc_trim = pd.read_excel(
        _projectDir + "/Output/2-TrimmoQC/QC_Table/" + _projectName + "_Trimmed_FastQC_Table.xlsx",
        index_col=0, header=[1, 2])

    # _currBatch_qc = join_levels(_currBatch_qc)
    # _currBatch_qc_trim = join_levels_mat3(_currBatch_qc_trim)
    _currBatch_qc_trim = join_levels(_currBatch_qc_trim)

    # print(_currBatch_qc_trim.columns)

    _currBatch_qc_trim = _currBatch_qc_trim.loc[:,
                         ['Sample', 'Overrepresented Sequences [% of Total Sequences]', '% of Adapter Content']]
    _currBatch_qc_trim.columns = ['Sample', 'Percent_Overrepresented_Seq_Trimmed', 'Percent_Adapter_Content_Trimmed']

    _currBatch_qc_trim.loc[:, 'Percent_Overrepresented_Seq_Trimmed'] = pd.to_numeric(
        _currBatch_qc_trim['Percent_Overrepresented_Seq_Trimmed'].str.strip("%"))
    _currBatch_qc_trim.loc[:, 'Percent_Adapter_Content_Trimmed'] = pd.to_numeric(
        _currBatch_qc_trim['Percent_Adapter_Content_Trimmed'].str.strip("%"))

    # print(_currBatch_qc_trim)

    _currBatch_qc = pd.merge(_currBatch_qc_untrim, _currBatch_qc_trim, on='Sample', how='left')
    # print(_currBatch_qc.columns)

    _currBatch_df = pd.merge(_currBatch_stat, _currBatch_qc, on='Sample', how='left')

    print(_currBatch_df)
    # print(_currBatch_df.iloc[0, 8])
    # print(_currBatch_df.columns)
    # print(_currBatch_df.Input_Size.mean())

    '''END of Current/Batch Input from USER}'''



    ## Adding last row in the dataframe with current library's mean values to be plotted on the final summary page
    # _currBatch_df.loc[:, "Percent_PostTrim"] = _retro_df.loc[:, "Percent_PostTrim"].clip(lower=60)
    _currBatch_summary_df = ["Library_Mean",
                             _currBatch_df.Input_Size.mean(),
                             _currBatch_df.Percent_PostTrim.mean(),
                             _currBatch_df.Num_Uniquely_Aligned.mean(),
                             _currBatch_df.Percent_Uniquely_Aligned.mean(),
                             _currBatch_df.Percent_Exonic.mean(),
                             _currBatch_df.Num_Uniquely_Aligned_rRNA.mean(),
                             _currBatch_df.Num_expressed_genes.mean(),
                             _currBatch_df.iloc[0, 8],
                             _currBatch_df.iloc[0, 9],
                             _currBatch_df.Percent_Overrepresented_Seq_Untrimmed.mean(),
                             _currBatch_df.Percent_Adapter_Content_Untrimmed.mean(),
                             _currBatch_df.Percent_Overrepresented_Seq_Trimmed.mean(),
                             _currBatch_df.Percent_Adapter_Content_Trimmed.mean()]

    print(_currBatch_summary_df)
    _currBatch_df.loc[len(_currBatch_df)] = _currBatch_summary_df

    print(_currBatch_df)

    '''{ REMOVE GC and HISTOGRAM data input for the FIRST ITERATION'''

    ## Read Gene Coverage Data
    _gc_df = pd.read_csv(_GC_data, index_col="Xaxis")

    ## Adding the Library Mean column at the end of the GC dataframe
    _gc_df["Library_Mean"] = _gc_df[_gc_df.columns].mean(axis=1)
    print(_gc_df)

    ## Read Histogram Distribution data
    _negBin_df = pd.read_csv(_hist_data, index_col=False)

    _negBin_df["Library_Mean"] = _negBin_df.iloc[:, 1:].mean(axis=1)
    print(_negBin_df)
    # print(_negBin_df.iloc[: ,1:])

    # print(next(_currBatch_df.iterrows())[1]['Sample'])

    '''END of GC and HISTOGRAM data USER input}'''

    '''{ BEGIN HISTORICAL data input/formatting '''

    ########### Based on USER input, load/format historical data from SCRIPT dataset
    ### or USER data (later stage) ###############

    ### Load the Metrics from masterDB
    _conn = create_connection(_masterDB)
    _cur = _conn.cursor()

    _master_stat = pd.read_sql_query("Select * from MasterMetrics", con=_conn)
    _conn.close()

    print(_master_stat['Date'])

    # _master_stat.loc[:, 'Date'] = pd.to_datetime(_master_stat['Date']).dt.date
    # _master_stat.loc[:, 'Date'] = pd.to_datetime(_master_stat['Date'], format='%Y-%m-%d')
    _master_stat['Date'] = pd.to_datetime(_master_stat['Date']).dt.date.astype('O')
    _master_stat['Date'] = pd.to_datetime(_master_stat['Date']).dt.date.astype('O')

    print(_master_stat.loc[2, 'Date'])

    # _master_stat.loc[:, 'Date'] = _master_stat.loc[:, 'Date'].apply(lambda t: t.strftime('%Y-%m-%d'))
    _master_stat.loc[:, 'Input_Size'] = _master_stat.loc[:, 'Input_Size'].apply(lambda j: (j / 1000000))

    print(_master_stat.columns)

    print(_master_stat.Date.dtype)
    print(_currBatch_df.Date.dtype)

    _currBatch_df_violin = _currBatch_df.drop('Date', axis=1)
    _master_stat_violin = _master_stat.drop('Date', axis=1)

    print(_currBatch_df_violin.columns)
    print(_master_stat_violin.columns)

    '''END of Historical data input/formatting}'''

    '''{BEGIN plotting process for Stage 1 - Stage 2'''

    ###### Remove the GC/HISTOGRAM plots for Stage 1 - Stage2
    ###### Create 6 panel grid instead of 8 panels

    ### Begin Plotting process
    _pdfObj = PdfPages(_plotDir + "/" + _projectName + "_retroPlot.pdf")

    for _tuple in _currBatch_df.itertuples():
        print(_tuple)

        fig = plt.figure(frameon=False)

        # Plotting figure 1: Input Size
        fig = plotHist_ipSize(_tuple, _currBatch_df, _master_stat, 1, fig)

        # Plotting figure 2: Percentage of Reads after Trimming
        fig = plotHist(_tuple, _currBatch_df, _master_stat, "Percent_PostTrim", "Trimming", 2, fig)

        # Plotting figure 3: Percentage of Uniquely Aligned Reads
        fig = plotHist(_tuple, _currBatch_df, _master_stat, "Percent_Uniquely_Aligned", "Alignment", 3, fig)

        # Plotting figure 4: Percentage of Reads Mapped to Exons
        fig = plotHist(_tuple, _currBatch_df, _master_stat, "Percent_Exonic", "Gene Exon Mapping", 4, fig)

        # Plotting figure 5: Scatter Plot of Number of Ribosomal RNA reads per Uniquely Aligned Reads
        fig = plotScatter(_tuple, _currBatch_df, _master_stat, 5, fig)

        # Plotting figure 6: Violin Plot for Contamination - % Adapter Content and % Overrepresented Sequences
        fig = plotViolin_dualAxis(_tuple, _currBatch_df_violin, _master_stat_violin, 6, fig)

        # Plotting figure 7: Gene Coverage plot for 5'-->3' coverage
        fig = plotGC(_tuple, _gc_df, 7, "GeneBody Coverage", fig)

        # Plotting Negative Binomial Histogram
        fig = plotNegBin(_tuple, _negBin_df, _currBatch_df, 8, "Distribution of Gene Expression", fig)

        if _tuple[1] == "Library_Mean":
            fig.suptitle('Sample : ' + _tuple[1], x=0.01, y=0.99, fontsize=6,
                         horizontalalignment='left', verticalalignment='top', fontweight='book', style='italic')

        else:
            fig.suptitle('SampleID : ' + _tuple[1] + '\t\t PatientID : ' + patient_mapper(_projectDir, _tuple[1]),
                         x=0.01, y=0.99, fontsize=6,
                         horizontalalignment='left', verticalalignment='top', fontweight='book', style='italic')

        fig.text(x=0.99, y=0.99, s=int(_tuple[0]) + 1, ha='right', va='top', fontsize=4)

        plt.subplots_adjust(hspace=0.7, wspace=0.2)

        _pdfObj.savefig(fig)

    _pdfObj.close()

    '''END plotting process for the first 6 panels}'''


def retroPlotter_main(_input_file):
    ### Read input file and load USER data
    _user_data = pd.read_csv(_input_file, sep = "\t")

    print(_user_data)






    return None




'''
def patient_mapper(_project_dir, _sample_id):
    _pat_dict = {}
    os.chdir(_project_dir)

    for _f in glob.glob("*DirMap.csv"):

        with open(_f, "r") as _patMap_file:
            _readr = csv.reader(_patMap_file, delimiter=',')
            for _row in _readr:
                _pat_dict[_row[1]] = _row[0]

    # print(_pat_dict[_sample_id])

    # print (_pat_dict)

    if _sample_id == "Library_Mean":
        return "Library_Mean"
    else:
        return _pat_dict[_sample_id]

'''





########## Remove SLURM orchestration and allow USER to run on native command line
## Python/anaconda3.6
## Matplotlib/3.0.0
## Pandas : ????

### Function for calling the RetroPlotter from the pipeline orchestration
def runRetro(_project_dir, _project_name, _gc_data, _hist_data, _plotOutDir):

    os.chdir(_project_dir + "/Output/4-htseqCounts/All_Htseq_scripts")

    script_name = "retroPlot_Caller.sh"
    script_f = open(script_name, "w+")

    script = ""
    script += "#!/bin/bash\n"
    script += "#SBATCH -A b1042\n"
    script += "#SBATCH -p genomics\n"
    script += "#SBATCH -t 24:00:00\n"
    script += "#SBATCH --mail-type=FAIL\n"
    script += "#SBATCH --mail-user=" + _notification_email + "\n"
    script += "#SBATCH --output=%x.%j.out\n"
    script += "#SBATCH --job-name=retroPlot_Wrapper\n"
    script += "#SBATCH -N 1\n"
    script += "#SBATCH -n 6\n\n"

    for module in required_modules:
        script += "module load " + module + "\n"

    script += "source activate pyMat3\n\n"

    script += "cd " + _pipeline_dir + "\n\n"

    script += "python helper_retroPlotter_DEV.py -pdir " + _project_dir + " -pname " + _project_name + " -gc " + _gc_data + " -hist " + _hist_data + " -out " + _plotOutDir + "\n\n"

    script_f.write(script)
    script_f.close()

    time.sleep(5)

    subprocess.call(['sbatch', script_name])

    print("RetroPlot script submitted for processing!\n\n")

    while not os.path.exists(_plotOutDir + "/" + _project_name + "_retroPlot.pdf"):

        print("Waiting for retro-Plotter to finish plotting. . .  .    .     .      .         .\n\n")
        time.sleep(5)
        pass
    else:
        time.sleep(300)
        return True





'''
if __name__ == "__main__":

    ''''''{BEGIN Take command line arguments for all input variables from the USER''''''

    ########## DESIGN new command line arguments for Stage 1 - Stage 2 development

    ### Run RetroPlotter from the pipeline using commandline arguments (Original Pipeline Method)

    parser = argparse.ArgumentParser(description="RetroPlotter Argument Parser.")

    # parser.add_argument("-v", "--version", help="Display U19 RNA-seq pipeline version.\n-v, --version\n", action="store_true")

    parser.add_argument("-pdir", "--project-dir", type=str, required=True,
                        help="[REQUIRED] Provide the absolute path of the project directory.\n -pdir [~/PROJECT-PATH/],\t--project-dir [~/PROJECT-PATH/]\n")
    parser.add_argument("-pname", "--project-name", type=str, required=True,
                        help="[REQUIRED] Provide the name of your project.\n -pname [PROJECT-NAME],\t--project-name [PROJECT-NAME]\n")
    parser.add_argument("-gc", "--genecoverage-data", type=str, required=True,
                        help="[REQUIRED] Provide the name of the Gene Coverage Data file containing the GC data to be plotted\n -gc [GC-DATA],\t--genecoverage-data [GC-DATA]\n")
    parser.add_argument("-hist", "--histogram-data", required=True, type=str,
                        help="[REQUIRED] Provide the name of the Histogram Data file containing the gene expression histogram data to be plotted.\n-hist [HIST-DATA],\t--histogram-data [HIST-DATA]\n")
    parser.add_argument("-out", "--output-dir", required=True, type=str,
                        help="[REQUIRED] Provide the name of the Output Directory for the plot.\n -out [OUTPUT-DIRNAME],\t--output-dir [OUTPUT-DIRNAME]\n")

    args = parser.parse_args()

    _pd = args.project_dir
    _pname = args.project_name
    _gc = args.genecoverage_data
    _histo = args.histogram_data
    _outdir = args.output_dir

    retroPlotter_mat3(_pd, _pname, _gc, _histo, _outdir)

    print("Successfully created the historical data plot.")


    ''''''END command line args from USER}''''''

'''

    

if __name__ == "__main__":
    
    
    ### Run RetroParser to take input using commandline arguments (USER Input)
    parser = argparse.ArgumentParser(description="RetroPlotter Argument Parser")

    parser.add_argument("-ip", "--input-file", type=str, required=True,
                        help="[REQUIRED] Provide the absolute path of the USER input file.\n -ip [~/INPUT-FILE-PATH/],\t--input-file [~/INPUT-FILE-PATH/]\n")


    parser.add_argument("-bgd", "--background-data", type=str, required=True,
                        help="[REQUIRED] Use the SCRIPT historical background data.\n -bgd [PROJECT-NAME],\t--project-name [PROJECT-NAME]\n")
    

    '''
    parser.add_argument("-gc", "--genecoverage-data", type=str, required=True,
                        help="[REQUIRED] Provide the name of the Gene Coverage Data file containing the GC data to be plotted\n -gc [GC-DATA],\t--genecoverage-data [GC-DATA]\n")
    parser.add_argument("-hist", "--histogram-data", required=True, type=str,
                        help="[REQUIRED] Provide the name of the Histogram Data file containing the gene expression histogram data to be plotted.\n-hist [HIST-DATA],\t--histogram-data [HIST-DATA]\n")
    parser.add_argument("-out", "--output-dir", required=True, type=str,
                        help="[REQUIRED] Provide the name of the Output Directory for the plot.\n -out [OUTPUT-DIRNAME],\t--output-dir [OUTPUT-DIRNAME]\n")
    '''
    

    args = parser.parse_args()

    _ip_file = args.input_file


    retroPlotter_main(_ip_file)


    

    '''
    ## Batch11   
    
    retroPlotter_mat3("/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11",
                      "SCRIPT_RNAseq_Batch_11",
                      "/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11/Output/3-STAROutput_Bam/GeneBody_Coverage/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv",
                      "/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11/Output/4-htseqCounts/Histogram_Plot/CPM_SCRIPT_RNAseq_Batch_11_final_count_bincounts_0.5.csv",
                      "/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11/Output/3-STAROutput_Bam/stats")
    
    '''