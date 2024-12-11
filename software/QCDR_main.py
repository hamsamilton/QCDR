import cProfile
import pstats
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import os
import math as math
import argparse
matplotlib.use('PDF')
import pickle
import matplotlib.pyplot as plt
import sys
import json
import seaborn as sns
import time
import shutil
import glob
import csv
import subprocess
import numpy as np
import pandas as pd
from helper_retroFunctions import *
from Sam_PyUtils import *
from sklearn.linear_model import LinearRegression
from scipy.stats import norm
from scipy import stats
from sklearn.preprocessing import LabelEncoder
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import MinMaxScaler
import matplotlib.offsetbox
from statsmodels.stats.weightstats import ztest

'''RetroPlotter caller function for reading data and passing it to individual plotters. Add option/flag for including GC/Hist and create 6-panel or 8-panel grid based on the flag passed to plotter functions'''

def QCDR_main(qry_filename    = '',
              op_folder       = '',
              _bgd_file       = '',
              _gc_file        = None,
              _hist_file      = None,
              cutoff_filename = None,
              fail_alpha      = .05,
              warn_alpha      = .1):

    #Make directory to store files in
    os.makedirs(op_folder,
                exist_ok = True)

    # Read query file and load USER data
    _user_df = pd.read_csv(qry_filename)
    _user_df = input_adapter(_user_df).adapt_input().input_df

    ## Read Background file 
    _bgd_df = pd.read_csv(_bgd_file)
    _bgd_df  = input_adapter(_bgd_df).adapt_input().input_df

    # Make standard cutoffs for warn/fail
    _fail_cutoffs = gen_cutoffs(bgd_df = _bgd_df,
                                alph   = fail_alpha)
    _warn_cutoffs = gen_cutoffs(bgd_df = _bgd_df,
                                alph   = warn_alpha)

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
    _figinfo["_bin_num"]           = 40
    _figinfo["_tick_size"]         = 4
    _figinfo["_subplot_rows"]      = _subplot_rows
    _figinfo["warn_alpha"]         = warn_alpha
    _figinfo["fail_alpha"]         = fail_alpha
    _figinfo["_ip_filename"]       = qry_filename
    _figinfo["_op_filename"]       = op_folder
    _figinfo["_gc_file"]           = _gc_file
    _figinfo["_hist_file"]         = _hist_file
    _figinfo["_bgd_filename"]      = _bgd_file
    _figinfo["_cutoff_filename"]   = cutoff_filename

    if cutoff_filename is not None:
        _manual_cutoffs = pd.read_excel(cutoff_filename)

        manual_cutoff_adaptr = manual_cutoff_adapter(_manual_cutoffs)
        manual_cutoff_adaptr.adapt_input()

        _man_warn_cutoff_dict = manual_cutoff_adaptr.man_cutoff_df['Warn'].to_dict()
        _man_fail_cutoff_dict = manual_cutoff_adaptr.man_cutoff_df['Fail'].to_dict()

        # for cutoffs with unspecified values, replace with the automatically generated cutoffs
        repl_missing_values_indict(_man_warn_cutoff_dict,_warn_cutoffs)
        repl_missing_values_indict(_man_fail_cutoff_dict,_fail_cutoffs)

        _warn_cutoffs = _man_warn_cutoff_dict
        _fail_cutoffs = _man_fail_cutoff_dict
        _warn_cutoffs["_alpha"] = _figinfo["warn_alpha"]
        _fail_cutoffs["_alpha"] = _figinfo["fail_alpha"]

    # add cutoff info
    _figinfo["_fail_cutoffs"] = _fail_cutoffs
    _figinfo["_warn_cutoffs"] = _warn_cutoffs

    # Read Gene Coverage Data
    if _gc_file is not None:
        _gc_df = pd.read_csv(_gc_file).iloc[:,1:]

        GCDeviances = EstimateDeviances(data_df = _gc_df).sum(axis = 0)

        _fail_GC_cutoff = CalcBootstrapBound(vec         = GCDeviances,
                                             alpha       = 2*fail_alpha,
                                             upper_lower = 'upper')
        _warn_GC_cutoff = CalcBootstrapBound(vec         = GCDeviances,
                                             alpha       = 2*warn_alpha,
                                             upper_lower = 'upper')
        _figinfo['_fail_cutoffs']['GC_cutoff'] = _fail_GC_cutoff
        _figinfo['_warn_cutoffs']['GC_cutoff'] = _warn_GC_cutoff
        GC_KS_pvals = stats.norm.sf(stats.zscore(GCDeviances))
        print('the GC deviances are')

        # Convert GC into KS vals and calculate the distribution to get a pvalue
        _figinfo["_gbc_pvals"]  = GC_KS_pvals
        _user_df["_gbc_pvals"]  = GC_KS_pvals
        _user_df = _user_df.set_index('Sample')
        print('the userdf was',_user_df,_user_df.index)
        _user_df['GBC_KSstats'] = GCDeviances
        _user_df = _user_df.reset_index()
        print('the _user_df is',_user_df)

        _figinfo["_gbc_exists"] = True
    else:
        _figinfo["_gbc_exists"] = False

    # Read Histogram data
    if _hist_file is not None:
        _negBin_df = CountsMatrixToGeneHist(df       = pd.read_excel(_hist_file),
                                            binsize  = .25,
                                            maxdepth = 18.5)
        # Preprocess raw counts table
        _data_df = _negBin_df.drop(['Bins'],
                                   axis = 1)
        _sum_df = _data_df.sum().round()
        _fail_numGene_cutoff = CalcBootstrapBound(vec         = _sum_df,
                                                  alpha       = 2*fail_alpha,
                                                  upper_lower = "lower")
        _warn_numGene_cutoff = CalcBootstrapBound(vec         = _sum_df,
                                                  alpha       = 2*warn_alpha,
                                                  upper_lower = "lower")
        _figinfo["_fail_cutoffs"]["_numGene_cutoff"] = _fail_numGene_cutoff
        _figinfo["_warn_cutoffs"]["_numGene_cutoff"] = _warn_numGene_cutoff
        _user_df["_hist_pvals"]  = calcHistPval(_data_df)
        _user_df["NumGenes"]  =  _data_df.sum().round().values
        _figinfo["_hist_pvals"]  = calcHistPval(_data_df)
        _figinfo["_hist_exists"] = True
    else:
        _figinfo["_fail_cutoffs"]["_numGene_cutoff"] = "None"
        _figinfo["_warn_cutoffs"]["_numGene_cutoff"] = "None"
        _figinfo["_hist_pvals"]  = None
        _figinfo["_hist_exists"] = False


    # Save the user_df
    _user_df.to_csv(op_folder + '/QCDR_ReportInfo.csv',
                    index = False)

    ###### Begin Plotting process ######
    pdf_output = op_folder + '/QCDR_Output.pdf'
    # Open the given PDF output file
    _pdfObj = PdfPages(pdf_output)

    # Create title page
    _title_fig = mkTitlePage(_figinfo)
    _pdfObj.savefig(_title_fig)
    plt.close()

    # how many tables do we need? 
    _summary_heatmap_data = mkQC_heatmap_data(_user_df,_figinfo)
    _summary_heatmap_data = pd.DataFrame(_summary_heatmap_data)
    # Save heatmap data
    _summary_heatmap_data["Sample"] = _user_df.Sample
    _summary_heatmap_data.to_csv(op_folder + '/Heatmapinfo.csv',
                                 index = False)
    my_range = list(range(0,len(_user_df),20))
    my_range.append(len(_user_df))

    for i in range(0,len(my_range)-1):
        new_rng = list(range(my_range[i],my_range[i+1]))
        _sub_df = _summary_heatmap_data.iloc[new_rng]
        _summary_heatmap_fig = mkQC_heatmap(_sub_df)

        _pdfObj.savefig(_summary_heatmap_fig)
        plt.close(_summary_heatmap_fig)

    for SampleName in _user_df["Sample"]:

        # Create empty figure
        fig = plt.figure(frameon=False)

        # Plotting figure 1: Input Size
        InputSize = ReadDepthHistPlotter(SampleName,_user_df,_bgd_df,1,_figinfo,fig)
        fig = InputSize.Figure

        # Plotting figure 2: Percentage of Reads after Trimming
        TrimmingPercent = TrimmingPlotter(SampleName, _user_df, _bgd_df,2,_figinfo,fig)
        fig = TrimmingPercent.Figure

        # Plotting figure 3: Percentage of Uniquely Aligned Reads
        Alignment = AlignmentPlotter(SampleName, _user_df, _bgd_df,3,_figinfo,fig)
        fig = Alignment.Figure

        # Plotting figure 4: Percentage of Reads Mapped to Exons
        ExonMapping = ExonMappingPlotter(SampleName, _user_df, _bgd_df,4,_figinfo,fig)
        fig = ExonMapping.Figure

        # Plotting figure 5: Scatter Plot of Number of Ribosomal RNA reads per Uniquely Aligned Reads
        fig = helper_retroFunctions.plotScatter_rRNA(SampleName, _user_df, _bgd_df, 5,_figinfo,fig)

        # Plotting figure 6: Violin Plot for Contamination - % Adapter Content and % Overrepresented Sequences
        fig = helper_retroFunctions.plotViolin_dualAxis(SampleName, _user_df, _bgd_df, 6,_figinfo,fig)

        # Plotting figure 7: Expression Distribution Plot
        if _hist_file is not None:
            fig = helper_retroFunctions.plotNegBin(SampleName,_user_df,_negBin_df,7,_figinfo,fig)

        # Plotting figure 8: Gene Body Coverage Plot
        if _gc_file is not None:
            fig = helper_retroFunctions.plotGC(SampleName,_user_df, _gc_df, 8,_figinfo,fig)

        # Add sample info at the top-left corner of the page
        fig.text(s                   = 'Sample : ' + SampleName,
                 x                   = 0.01,
                 y                   = 0.99,
                 fontsize            = 6,
                 horizontalalignment = 'left',
                 verticalalignment   = 'top',
                 fontweight          = 'book',
                 style               = 'italic')
        fig.text(s                   = "Batch : " + _user_df.loc[_user_df['Sample'] == SampleName,'Batch'].iloc[0],
                 x                   = 0.99,
                 y                   = 0.99,
                 fontsize            = 6,
                 horizontalalignment = 'right',
                 verticalalignment   = 'top',
                 fontweight          = 'book',
                 style               = 'italic')

        plt.subplots_adjust(left   = .07,
                            right  = .93,
                            bottom = .05,
                            top    = .9,
                            hspace = .72,
                            wspace = 0.25)

        _pdfObj.savefig(fig)
        plt.close(fig)
    _pdfObj.close()

    return None


if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    ### Run RetroParser to take input using commandline arguments (USER Input)
    parser = argparse.ArgumentParser(description = "QCDR Argument Parser")

    parser.add_argument("-qry",
                        "--query-filename",
                        type     = os.path.abspath,
                        required = True,
                        help     = """[REQUIRED] Provide the location of the USER input file.
                                      -qry [~/QUERY-FILE-PATH/]""")

    parser.add_argument("-out",
                        "--output-folder",
                        required = False,
                        type     = os.path.abspath,
                        default  = "QCDR_Outputs",
                        help     = """[OPTIONAL] Where to save outputs.
                                      -out [OUTPUT-FILENAME.pdf],\t--output-dir [OUTPUT-FOLDER]""")

    parser.add_argument("-bgd",
                        "--background-data",
                        required = False,
                        help     ="[OPTIONAL] Location of background data to contextualize query dataset in")

    parser.add_argument("-gc",
                        "--genecoverage-data",
                        type     = str,
                        required = False,
                        help     = """[OPTIONAL] Provide the name of the Gene Coverage Data file
                                      -gc [GC-DATA],\t--genecoverage-data [GC-DATA]""")

    parser.add_argument("-hist",
                        "--histogram-data",
                        required = False,
                        type     = str,
                        help     = """[OPTIONAL] Provide a raw count matrix to generate a read count histogram
                                       -hist [HIST-DATA],\t--histogram-data [HIST-DATA]""")

    parser.add_argument("-ctf",
                        "--cutoff-data",
                        required = False,
                        default  = None,
                        help     = """[OPTIONAL] Location of table with cutoffs for warn fail cutoffs.
                                      -ctf [CUTOFF_PATH],\t --cutoff [CUTOFF_PATH]""")

    parser.add_argument("-fla",
                        "--failalpha",
                        required = False,
                        default  = .05,
                        help     = """[OPTIONAL] Provide an alpha cutoff for determining FAIL flags
                                      -fla [FAIL_ALPHA] \t --failalpha [FAIL_ALPHA]""")

    parser.add_argument("-wrna",
                        "--warnalpha",
                        required = False,
                        type     = float,
                        default  = .1,
                        help     ="""[OPTIONAL] Provide an alpha cutoff for warn.
                                     -wrna [WARN_ALPHA],\t--warnalpha [WARN_ALPHA]""")

    args = parser.parse_args()

    qry_filename    = args.query_filename
    op_folder       = args.output_folder
    gc_file         = args.genecoverage_data
    hist_file       = args.histogram_data
    bgd_filename    = args.background_data if args.background_data else qry_filename
    cutoff_filename = args.cutoff_data
    fail_alpha      = float(args.failalpha)
    warn_alpha      = float(args.warnalpha)

    print(f"Query Data : {qry_filename}")
    print(f"Output File : {op_folder}")
    print(f"Background Data : {bgd_filename}")
    print(f"cutoffs_provided : {cutoff_filename}")
    print(f"failalpha : {fail_alpha}")
    print(f"warnalpha : {warn_alpha}")

    QCDR_main(qry_filename    = qry_filename,
              op_folder       = op_folder,
              _bgd_file       = bgd_filename,
              _gc_file        = gc_file,
              _hist_file      = hist_file,
              cutoff_filename = cutoff_filename,
              fail_alpha      = fail_alpha,
              warn_alpha      = warn_alpha)

    profiler.disable()

    ps = pstats.Stats(profiler,
                      stream = sys.stdout).sort_stats('cumulative')
    ps.print_stats(20)
    print('Done')
