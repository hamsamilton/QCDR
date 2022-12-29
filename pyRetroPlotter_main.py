import matplotlib
import matplotlib.gridspec as gridspec
import os
import argparse
matplotlib.use('PDF')
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

_script_bgdFile = "/projects/b1063/Gaurav/pyRetroPlotter/data/SCRIPTretro_masterStatistics_allBatches.csv"
_masterDB = "/projects/b1079/tools/U19_masterDB_dev/U19_masterDB.sqlite"
_pipeline_dir = "/projects/b1079/U19_RNAseq_Pipeline_SLURM_v01"
required_modules = ['python/anaconda']
_notification_email = "samuelhamilton2024@u.northwestern.edu"

# Replaces missing values of a dictionary with corresponding values of a second dictionary with matching keys
def repl_missing_values_indict(_indict,_repldict):
    for key, value in _indict.items():
    # If the value is NaN, replace it with the corresponding automatically generated keys
        if value != value:
            _indict.update({key: _repldict[key]})

'''RetroPlotter caller function for reading data and passing it to individual plotters. Add option/flag for including GC/Hist and create 6-panel or 8-panel grid based on the flag passed to plotter functions'''

def retroPlotter_main(_input_file, _output_file, _bgd_file, _gc_file,_hist_file):

    print(type(_fail_alpha))
    # Read input file and load USER data
    _user_df = pd.read_csv(_input_file, sep=",")

    # Add last row in the USER df with current batch's mean values for the final summary page
    _batch_summary_df = ["Batch_Mean",
                         _user_df.Input_Size.mean(),
                         _user_df.Percent_PostTrim.mean(),
                         _user_df.Num_Uniquely_Aligned.mean(),
                         _user_df.Percent_Uniquely_Aligned.mean(),
                         _user_df.Percent_Exonic.mean(),
                         _user_df.Num_Uniquely_Aligned_rRNA.mean(),
                         _user_df.Percent_Overrepresented_Seq_Untrimmed.mean(),
                         _user_df.Percent_Adapter_Content_Untrimmed.mean(),
                         _user_df.Percent_Overrepresented_Seq_Trimmed.mean(),
                         _user_df.Percent_Adapter_Content_Trimmed.mean(),
                         _user_df.iloc[0, 11],
                         _user_df.iloc[0, 12]]

    ## Add Batch mean as the last row of the USER dataframe
    _user_df.loc[len(_user_df)] = _batch_summary_df

    ## Read Background file and load HISTORICAL background data
    _bgd_df = pd.read_csv(_bgd_file, sep=",")

    # Convert Input_Size to per Million (/1000000) for ease of plotting first panel
    _bgd_df.loc[:, 'Input_Size'] = _bgd_df.loc[:, 'Input_Size'].apply(lambda j: (j / 1000000))
    # Make standard cutoffs for warn/fail
    _fail_cutoffs = helper_retroFunctions.gen_cutoffs(_bgd_df = _bgd_df,_alph = _fail_alpha)

    _warn_cutoffs = helper_retroFunctions.gen_cutoffs(_bgd_df = _bgd_df,_alph = _warn_alpha)
    
    if _cutoff_filename != False:
        _manual_cutoffs = pd.read_excel(_cutoff_filename)
        print(_manual_cutoffs.to_string())
        _man_warn_cutoff_dict = _manual_cutoffs.set_index('cutoff')['Warn'].to_dict()
        _man_fail_cutoff_dict = _manual_cutoffs.set_index('cutoff')['Fail'].to_dict()
 
        # for cutoffs with unspecified values, replace with the automatically generated cutoffs

        repl_missing_values_indict(_man_warn_cutoff_dict,_warn_cutoffs)
        repl_missing_values_indict(_man_fail_cutoff_dict,_fail_cutoffs)

        _warn_cutoffs = _man_warn_cutoff_dict
        _fail_cutoffs = _man_fail_cutoff_dict
    
    ## Read Gene Coverage Data
    _gc_df = pd.read_csv(_gc_file, index_col="Xaxis")

    ## Adding the Library Mean column at the end of the GC dataframe
    _gc_df["Batch_Mean"] = _gc_df[_gc_df.columns].mean(axis=1)

    # Read Histogram data
    _negBin_df = pd.read_csv(_hist_file, index_col=False)
    _negBin_df["Batch_Mean"] = _negBin_df.iloc[:, 1:].mean(axis=1)

    ###### Begin Plotting process ######

    ## Open the given PDF output file
    _pdfObj = PdfPages(_output_file)

    for _tuple in _user_df.itertuples():
        print("THIS IS THE IN TUPLE",_tuple)

        # Create empty figure
        fig = plt.figure(frameon=False)

        # Plotting figure 1: Input Size
        fig = helper_retroFunctions.plotHist_ipSize(_tuple, _user_df, _bgd_df, 1,_fail_cutoffs["_ipReads_cutoff"],_warn_cutoffs["_ipReads_cutoff"],fig)

        # Plotting figure 2: Percentage of Reads after Trimming
        fig = helper_retroFunctions.plotHist_trimming(_tuple, _user_df, _bgd_df, "Percent_PostTrim", "Trimming", 2,_fail_cutoffs["_trimmedReads_cutoff"],_warn_cutoffs["_trimmedReads_cutoff"],fig)

        # Plotting figure 3: Percentage of Uniquely Aligned Reads
        fig = helper_retroFunctions.plotHist_alignment(_tuple, _user_df, _bgd_df, "Percent_Uniquely_Aligned", "Alignment", 3,_fail_cutoffs["_uniqAligned_cutoff"],_warn_cutoffs["_uniqAligned_cutoff"],fig)

        # Plotting figure 4: Percentage of Reads Mapped to Exons
        fig = helper_retroFunctions.plotHist_exonMapping(_tuple, _user_df, _bgd_df, "Percent_Exonic", "Gene Exon Mapping", 4,_fail_cutoffs["_exonMapping_cutoff"],_warn_cutoffs["_exonMapping_cutoff"], fig)

        # Plotting figure 5: Scatter Plot of Number of Ribosomal RNA reads per Uniquely Aligned Reads
        fig = helper_retroFunctions.plotScatter_rRNA(_tuple, _user_df, _bgd_df, 5,_fail_cutoffs["_riboScatter_cutoff"],_warn_cutoffs["_riboScatter_cutoff"],fig)

        # Plotting figure 6: Violin Plot for Contamination - % Adapter Content and % Overrepresented Sequences
        fig = helper_retroFunctions.plotViolin_dualAxis(_tuple, _user_df, _bgd_df, 6,_fail_cutoffs["_violin_cutoff_overrep_untrimmed"],_fail_cutoffs["_violin_cutoff_adapter_untrimmed"],
        _warn_cutoffs["_violin_cutoff_overrep_untrimmed"],_warn_cutoffs["_violin_cutoff_adapter_untrimmed"],_fail_cutoffs["_violin_cutoff_overrep_trimmed"],_fail_cutoffs["_violin_cutoff_adapter_trimmed"],
        _warn_cutoffs["_violin_cutoff_overrep_trimmed"],_warn_cutoffs["_violin_cutoff_adapter_trimmed"],fig)

        # Plotting figure 7: GeneBody Coverage Distribution Plot
        fig = helper_retroFunctions.plotGC(_tuple, _gc_df, 7, "GeneBody Coverage Distribution",_fail_alpha,_warn_alpha,fig)

        # Plotting figure 8: Gene Expression Distribution Plot
        fig = helper_retroFunctions.plotNegBin(_tuple,_negBin_df,_user_df,8,"Distribution of Gene Expression",_fail_alpha,_warn_alpha,fig)

        # Add sample name at the top-left corner of the page
        fig.suptitle('Sample : ' + _tuple[1] + " Batch : " + _tuple[13], x=0.01, y=0.99, fontsize=6,
                     horizontalalignment='left', verticalalignment='top', fontweight='book',style = 'italic')
        fig.text(s= ('Warn cutoffs: Default alpha =' + str(_warn_alpha) +
                    '| Sequencing Depth = ' + str(_warn_cutoffs["_ipReads_cutoff"]) +
                    '| Trimming = ' + str(_warn_cutoffs["_trimmedReads_cutoff"]) +  
                    '| Alignment = ' + str(_warn_cutoffs["_uniqAligned_cutoff"]) +  
                    '| Gene Exon Mapping = ' + str(_warn_cutoffs["_exonMapping_cutoff"]) + 
                    '| Ribosomal RNA = ' + str(round(_warn_cutoffs["_riboScatter_cutoff"],3)) + 
                    '| Sequence Contamination = ' + str(_warn_cutoffs["_violin_cutoff_adapter_untrimmed"]) + 
                    '| Gene Body Coverage = ' + str(_warn_alpha) + 
                    '| Distribution of Gene Expression = ' + str(_warn_alpha)) , 
                    x = .01, y = .965, fontsize = 3,
                     ha='left', va='top',fontweight='book', style = 'italic')
        fig.text(s= ('Fail cutoffs: Default alpha =' + str(_fail_alpha) +
                    '| Sequencing Depth = ' + str(_fail_cutoffs["_ipReads_cutoff"]) +
                    '| Trimming = ' + str(_fail_cutoffs["_trimmedReads_cutoff"]) +  
                    '| Alignment = ' + str(_fail_cutoffs["_uniqAligned_cutoff"]) +  
                    '| Gene Exon Mapping = ' + str(_fail_cutoffs["_exonMapping_cutoff"]) + 
                    '| Ribosomal RNA = ' + str(round(_fail_cutoffs["_riboScatter_cutoff"],3)) + 
                    '| Sequence Contamination = ' + str(_fail_cutoffs["_violin_cutoff_adapter_untrimmed"]) + 
                    '| Gene Body Coverage = ' + str(_fail_alpha) + 
                    '| Distribution of Gene Expression = ' + str(_fail_alpha)) , 
                    x = .01, y = .95, fontsize = 3,
                     ha='left', va='top',fontweight='book', style = 'italic')
        # Add page number at the top-right corner of the page
        fig.text(x=0.99, y=0.99, s=int(_tuple[0]) + 1, ha='right', va='top', fontsize=4)
            
        plt.subplots_adjust(hspace=0.7, wspace=0.2)

        _pdfObj.savefig(fig)
        plt.close(fig) 
    _pdfObj.close()

    return None


if __name__ == "__main__":
    
    ### Run RetroParser to take input using commandline arguments (USER Input)
    parser = argparse.ArgumentParser(description="RetroPlotter Argument Parser")

    parser.add_argument("-ip", "--input-filename", type=os.path.abspath, required=True,
                        help="[REQUIRED] Provide the name of the USER input file.\n -ip [~/INPUT-FILE-PATH/],\t--input-file [~/INPUT-FILE-PATH/]\n")

    parser.add_argument("-out", "--output-filename", required=False, type=os.path.abspath, default = "pyRetroPlotter_OutputPlot.pdf",
                        help="[REQUIRED] Provide the name of the PDF output file for the plot.\n -out [OUTPUT-FILENAME.pdf],\t--output-dir [OUTPUT-FILENAME.pdf]\n")

    parser.add_argument("-bgd", "--background-data", required=False, default = "/projects/b1063/Gaurav/pyRetroPlotter/data/SCRIPTretro_masterStatistics_allBatches.csv", 
                        help="[OPTIONAL] Use the SCRIPT historical background data or specify your own.")

    parser.add_argument("-gc", "--genecoverage-data", type=str, required=False,
                        help="[OPTIONAL] Provide the name of the Gene Coverage Data file containing the GC data to be plotted\n -gc [GC-DATA],\t--genecoverage-data [GC-DATA]\n")

    parser.add_argument("-hist", "--histogram-data", required=False, type=str,
                        help="[OPTIONAL] Provide the name of the Histogram Data file containing the gene expression histogram data to be plotted.\n-hist [HIST-DATA],\t--histogram-data [HIST-DATA]\n")
    
    parser.add_argument("-ctf", "--cutoffs", required=False, default = False,help="[OPTIONAL] Provide optional cutoffs for warn fail cutoffs.\n -ctf [CUTOFF_PATH],\t --cutoff [CUTOFF_PATH] \n")

    parser.add_argument("-fla", "--failalpha", required=False, default = .9,help="[OPTIONAL] Provide an alpha cutoff for failure.\n -fla [FAIL_ALPHA],\t--fail-alpha [FAIL_ALPHA]\n") 
    
    parser.add_argument("-wrna", "--warnalpha", required=False,type=float, default = .8,help="[OPTIONAL] Provide an alpha cutoff for warn.\n -wrna [WARN_ALPHA],\t--warnalpha [WARN_ALPHA]\n")
 
    args = parser.parse_args()

    _ip_filename     = args.input_filename
    _op_filename     = args.output_filename
    _gc_file         = args.genecoverage_data
    _hist_file       = args.histogram_data
    _bgd_filename    = args.background_data
    _cutoff_filename = args.cutoffs
    _fail_alpha      = float(args.failalpha)
    _warn_alpha      = float(args.warnalpha)    
   
    print(f"Input File : {_ip_filename}")
    print(f"Output File : {_op_filename}")
    print(f"Background File : {_bgd_filename}")
    print(f"cutoffs_provided : {_cutoff_filename}")
    print(f"failalpha : {_fail_alpha}")
    print(f"warnalpha : {_warn_alpha}")
    
    retroPlotter_main(_ip_filename, _op_filename, _bgd_filename,_gc_file,_hist_file)


    

    '''
    ## Batch11   
    
    retroPlotter_mat3("/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11",
                      "SCRIPT_RNAseq_Batch_11",
                      "/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11/Output/3-STAROutput_Bam/GeneBody_Coverage/SCRIPT_RNAseq_Batch_11_GeneCoverageData.csv",
                      "/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11/Output/4-htseqCounts/Histogram_Plot/CPM_SCRIPT_RNAseq_Batch_11_final_count_bincounts_0.5.csv",
                      "/projects/b1042/WinterLab/SCRIPT_ComplexityAnalysis/Datasets/SCRIPT_RNAseq_Batch_11/Output/3-STAROutput_Bam/stats")
    '''
