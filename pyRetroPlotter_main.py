import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import os
import math as math
import argparse
matplotlib.use('PDF')
import seaborn
import pickle
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
from helper_retroFunctions import *
from Sam_PyUtils import *
from sklearn.linear_model import LinearRegression
from itertools import zip_longest
from scipy.stats import norm
from scipy import stats
from sklearn.preprocessing import LabelEncoder
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import MinMaxScaler
import matplotlib.offsetbox
from datetime import datetime
# from pandas.plotting import register_matplotlib_converters
from statsmodels.stats.weightstats import ztest

'''RetroPlotter caller function for reading data and passing it to individual plotters. Add option/flag for including GC/Hist and create 6-panel or 8-panel grid based on the flag passed to plotter functions'''

def retroPlotter_main(_input_file, _output_file, _bgd_file, _gc_file,_hist_file):
    
     # Read input file and load USER data
    _user_df = pd.read_csv(_input_file, sep=",")
 
    input_adaptr=input_adapter(_user_df)   
    input_adaptr.adapt_input()
    _user_df = input_adaptr.input_df
    print("this is the user_df",_user_df)
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

    input_adaptr=input_adapter(_bgd_df)   
    input_adaptr.adapt_input()
    _bgd_df = input_adaptr.input_df
    print("this is the bgd_df",_bgd_df)


    # Convert Input_Size to per Million (/1000000) for ease of plotting first panel
    _bgd_df.loc[:, 'Input_Size'] =  _bgd_df.loc[:, 'Input_Size'].apply(lambda j: (j / 1000000))
    _user_df.loc[:, 'Input_Size'] = _user_df.loc[:, 'Input_Size'].apply(lambda j: (j / 1000000))

    # Make standard cutoffs for warn/fail
    _fail_cutoffs = gen_cutoffs(bgd_df = _bgd_df,alph = _fail_alpha)
    _warn_cutoffs = gen_cutoffs(bgd_df = _bgd_df,alph = _warn_alpha)
    
    # add an additional row if the gc or hist data was added
    if _gc_file is None and _hist_file is None:
        _subplot_rows = 3
    else:
        _subplot_rows = 4
    
    # create dictionarys to store values in 2 pass 2 helper functions
    _figinfo = {}
    _figinfo["_fail_color"]        = "red"
    _figinfo["_warn_color"]        = "goldenrod"
    _figinfo["_curr_sample_color"] = "lightseagreen"
    _figinfo["_title_size"]        = 6
    _figinfo["_label_size"]        = 5
    _figinfo["_legend_size"]       = 3
    _figinfo["_tick_size"]         = 4
    _figinfo["_subplot_rows"]      = _subplot_rows
    _figinfo["_warn_alpha"]        = _warn_alpha
    _figinfo["_fail_alpha"]        = _fail_alpha
    _figinfo["_bin_num"]           = 40
    _figinfo["_ip_filename"]       = _ip_filename
    _figinfo["_op_filename"]       = _op_filename
    _figinfo["_gc_file"]           = _gc_file
    _figinfo["_hist_file"]         = _hist_file
    _figinfo["_bgd_filename"]      = _bgd_filename
    _figinfo["_cutoff_filename"]   = _cutoff_filename

    if _cutoff_filename != False:
        _manual_cutoffs = pd.read_excel(_cutoff_filename)
        
        manual_cutoff_adaptr = manual_cutoff_adapter(_manual_cutoffs)
        manual_cutoff_adaptr.adapt_input()
        
        _man_warn_cutoff_dict = manual_cutoff_adaptr.man_cutoff_df['Warn'].to_dict()
        _man_fail_cutoff_dict = manual_cutoff_adaptr.man_cutoff_df['Fail'].to_dict()
 
        # for cutoffs with unspecified values, replace with the automatically generated cutoffs

        repl_missing_values_indict(_man_warn_cutoff_dict,_warn_cutoffs)
        repl_missing_values_indict(_man_fail_cutoff_dict,_fail_cutoffs)

        _warn_cutoffs = _man_warn_cutoff_dict
        _fail_cutoffs = _man_fail_cutoff_dict
        _warn_cutoffs["_alpha"] = _figinfo["_warn_alpha"]
        _fail_cutoffs["_alpha"] = _figinfo["_fail_alpha"]
    # add cutoff info
    
    _figinfo["_fail_cutoffs"] = _fail_cutoffs
    _figinfo["_warn_cutoffs"] = _warn_cutoffs
    
    # Read Gene Coverage Data
    if _gc_file is not None:
        _gc_df = pd.read_csv(_gc_file, index_col="Xaxis")
        
        # Adding the Library Mean column at the end of the GC dataframe
        _gc_df["Batch_Mean"] = _gc_df[_gc_df.columns].mean(axis=1)
        
        gc_KSvals = GC_KSstats(_gc_df)
        
        # Calculate pvalues for plotting by the summary heatmap
        _gc_vals = stats.norm.sf(stats.zscore(gc_KSvals))

        # Convert the distribution of stuff
        _figinfo["_gbc_pvals"] = _gc_vals
        _user_df["_gbc_pvals"]  = _gc_vals
        _figinfo["_gbc_exists"] = True
    else:
        _figinfo["_gbc_exists"] = False

    # Read Histogram data
    if _hist_file is not None:
        _negBin_df = pd.read_csv(_hist_file, index_col=False)
        _negBin_df["Batch_Mean"] = _negBin_df.iloc[:, 1:].mean(axis=1)

        _data_df = _negBin_df.drop(['Unnamed: 0'], axis=1)
        _sum_df = _data_df.sum().round()
        _fail_numGene_cutoff = get_ci_bound(vec = _sum_df,
                                                alpha = 2*_fail_alpha,
                                            upper_lower = "lower")
        _warn_numGene_cutoff = get_ci_bound(vec = _sum_df,
                                                alpha = 2*_warn_alpha,
                                            upper_lower = "lower")
        _figinfo["_fail_cutoffs"]["_numGene_cutoff"] = '{:.0f}'.format(_fail_numGene_cutoff) 
        _figinfo["_warn_cutoffs"]["_numGene_cutoff"] = '{:.0f}'.format(_warn_numGene_cutoff)
        _user_df["_hist_pvals"]  = calcHistPval(_negBin_df) 
        _figinfo["_hist_pvals"]  = calcHistPval(_negBin_df)
        _figinfo["_hist_exists"] = True
    else:
        _figinfo["_fail_cutoffs"]["_numGene_cutoff"] = "None"
        _figinfo["_warn_cutoffs"]["_numGene_cutoff"] = "None"
        _figinfo["_hist_pvals"]  = None 
        _figinfo["_hist_exists"] = False
    ###### Begin Plotting process ######

    # Open the given PDF output file
    _pdfObj = PdfPages(_output_file)
   
    # Create title page
    _title_fig = mkTitlePage(_figinfo) 
    _pdfObj.savefig(_title_fig)
    plt.close()

    # how many tables do we need? 
    _summary_heatmap_data = mkQC_heatmap_data(_user_df,_figinfo)
    _summary_heatmap_data = pd.DataFrame(_summary_heatmap_data)
    _summary_heatmap_data["Sample"] = _user_df.Sample
    my_range = list(range(0,len(_user_df),20))
    my_range.append(len(_user_df))
    for i in range(0,len(my_range)-1):
        new_rng = list(range(my_range[i],my_range[i+1]))
        _sub_df = _summary_heatmap_data.iloc[new_rng]
        _summary_heatmap_fig = mkQC_heatmap(_sub_df)
    
        _pdfObj.savefig(_summary_heatmap_fig)
        plt.close(_summary_heatmap_fig)   
 
    for _tuple in _user_df.itertuples():

        # Create empty figure
        fig = plt.figure(frameon=False)

        # Plotting figure 1: Input Size
        fig = helper_retroFunctions.plotHist_ipSize(_tuple, _user_df, _bgd_df, 1,_figinfo,fig)

        # Plotting figure 2: Percentage of Reads after Trimming
        fig = helper_retroFunctions.plotHist_trimming(_tuple, _user_df, _bgd_df, "Percent_PostTrim", "Trimming", 2,_figinfo,fig)

        # Plotting figure 3: Percentage of Uniquely Aligned Reads
        fig = helper_retroFunctions.plotHist_alignment(_tuple, _user_df, _bgd_df, "Percent_Uniquely_Aligned", "Alignment", 3,
                                                       _figinfo,fig)

        # Plotting figure 4: Percentage of Reads Mapped to Exons
        fig = helper_retroFunctions.plotHist_exonMapping(_tuple, _user_df, _bgd_df, "Percent_Exonic", "Exon Mapping", 4,
                                                         _figinfo, fig)

        # Plotting figure 5: Scatter Plot of Number of Ribosomal RNA reads per Uniquely Aligned Reads
        fig = helper_retroFunctions.plotScatter_rRNA(_tuple, _user_df, _bgd_df, 5,_figinfo,fig)
 
        # Plotting figure 6: Violin Plot for Contamination - % Adapter Content and % Overrepresented Sequences
        fig = helper_retroFunctions.plotViolin_dualAxis(_tuple, _user_df, _bgd_df, 6,_figinfo,fig)

        # Plotting figure 7: Expression Distribution Plot
        if _gc_file is not None:
            fig = helper_retroFunctions.plotNegBin(_tuple,_negBin_df,_user_df,7,"Gene Expression",_figinfo,fig)           

        # Plotting figure 8: Gene Body Coverage Plot
        if _hist_file is not None:
            fig = helper_retroFunctions.plotGC(_tuple, _gc_df, 8, "GeneBody Coverage",_figinfo,fig)

        # Add sample info at the top-left corner of the page
        fig.text(s='Sample : ' + _tuple[1], x=0.01, y=0.99, fontsize=6,
                     horizontalalignment='left', verticalalignment='top', fontweight='book',style = 'italic')
        fig.text(s ="Batch : " + _tuple[13], x=0.99, y=0.99, fontsize=6,
                     horizontalalignment='right', verticalalignment='top', fontweight='book',style = 'italic')
            
        plt.subplots_adjust(left = .05,right = .95, bottom = .05, top = .9,hspace=.72, wspace=0.25)

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

    parser.add_argument("-fla", "--failalpha", required=False, default = .05,help="[OPTIONAL] Provide an alpha cutoff for failure.\n -fla [FAIL_ALPHA],\t--fail-alpha [FAIL_ALPHA]\n") 
    
    parser.add_argument("-wrna", "--warnalpha", required=False,type=float, default = .1,help="[OPTIONAL] Provide an alpha cutoff for warn.\n -wrna [WARN_ALPHA],\t--warnalpha [WARN_ALPHA]\n")
 
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
