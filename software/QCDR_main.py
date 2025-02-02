import cProfile
import pstats
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import os
import argparse
matplotlib.use('PDF')
import pickle
import matplotlib.pyplot as plt
import sys
import json
import time
import shutil
import glob
import csv
import subprocess
import numpy as np
import pandas as pd
from helper_retroFunctions import *
from scipy.stats import norm
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.offsetbox

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
    _user_df = read_file(qry_filename)
    _user_df = input_adapter(_user_df).adapt_input().input_df

    ## Read Background file 
    _bgd_df = read_file(_bgd_file)
    _bgd_df  = input_adapter(_bgd_df).adapt_input().input_df
    MetricInfo = InitMetricInfo()

    # Read Gene Coverage Data
    if _gc_file is not None:
        _gc_df = read_file(_gc_file).iloc[:,1:]

        GCDeviances = EstimateDeviances(data_df = _gc_df).sum(axis = 0)
        #GCDeviances = GC_KSstats(_coverage_df = _gc_df)

        GCDeviances = pd.Series(GCDeviances,
                                index = _gc_df.columns,
                                name  = 'GBC_KSstats')

        _user_df = _user_df.set_index('Sample')
        _bgd_df  = _bgd_df.set_index('Sample')

        _user_df["GBC_KSstats"] = _user_df.index.map(GCDeviances).fillna(0) # Not sure if this should be filling NAs
        _bgd_df["GBC_KSstats"] = _bgd_df.index.map(GCDeviances).fillna(0)

        _user_df = _user_df.reset_index()
        _bgd_df = _bgd_df.reset_index()

    # Read Histogram data
    if _hist_file is not None:
        _negBin_df = CountsMatrixToGeneHist(df       = read_file(_hist_file),
                                            binsize  = .25,
                                            maxdepth = 18.5)
        # Preprocess raw counts table
        _data_df = _negBin_df.drop(['Bins'],
                                   axis = 1)
        _sum_df = _data_df.sum().round()
        NumGenes_series = pd.Series(_sum_df.values,
                                    index = _data_df.columns,
                                    name  = "NumGenes")

        _user_df = _user_df.set_index("Sample")
        _bgd_df = _bgd_df.set_index("Sample")
        _user_df["NumGenes"] = _user_df.index.map(NumGenes_series).fillna(0)
        _bgd_df["NumGenes"]  = _bgd_df.index.map(NumGenes_series).fillna(0)

        _user_df = _user_df.reset_index()
        _bgd_df = _bgd_df.reset_index()

    # Make standard cutoffs for warn/fail
    _fail_cutoffs = CutoffCalculator(bgd_df      = _bgd_df,
                                     alph        = fail_alpha).gen_cutoffs()
    _warn_cutoffs = CutoffCalculator(bgd_df      = _bgd_df,
                                     alph        = warn_alpha).gen_cutoffs()

    # Add cutoffs to the cutoffinfo
    # Add _fail_cutoffs and _warn_cutoffs as new columns
    MetricInfo["Fail_Cutoff"] = MetricInfo["Metric"].map(_fail_cutoffs)
    MetricInfo["Warn_Cutoff"] = MetricInfo["Metric"].map(_warn_cutoffs)

    print(MetricInfo)

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
    _figinfo['_MetricInfo']        = MetricInfo

    if cutoff_filename is not None:
        _manual_cutoffs = read_file(cutoff_filename)

        manual_cutoff_adaptr = manual_cutoff_adapter(_manual_cutoffs)
        manual_cutoff_adaptr.adapt_input()

        _man_warn_cutoff_dict = manual_cutoff_adaptr.input_df['Warn'].to_dict()
        _man_fail_cutoff_dict = manual_cutoff_adaptr.input_df['Fail'].to_dict()

        # for cutoffs with unspecified values, replace with the automatically generated cutoffs
        repl_missing_values_indict(_man_warn_cutoff_dict,_warn_cutoffs)
        repl_missing_values_indict(_man_fail_cutoff_dict,_fail_cutoffs)

        MetricInfo["Fail_Cutoff"] = MetricInfo["Metric"].map(_man_fail_cutoff_dict)
        MetricInfo["Warn_Cutoff"] = MetricInfo["Metric"].map(_man_warn_cutoff_dict)

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
    _summary_heatmap_data.to_csv(op_folder + '/Heatmapinfo.csv',
                                 index = False)

    # Subset the heatmap data into smaller heatmaps for plotting
    my_range = list(range(0,len(_user_df),20))
    my_range.append(len(_user_df)) # how many rows will be used for plotting

    # Get subset and plot
    for i in range(0,len(my_range)-1):
        new_rng = list(range(my_range[i],my_range[i+1]))
        _sub_df = _summary_heatmap_data.iloc[new_rng]
        _summary_heatmap_fig = mkQC_heatmap(_sub_df)

        _pdfObj.savefig(_summary_heatmap_fig)
        plt.close(_summary_heatmap_fig)

    for SampleName in _user_df["Sample"]:
        print(SampleName)

        # Create empty figure
        fig = plt.figure(frameon=False)

        # Plotting figure 1: Input Size
        if pd.notna(MetricInfo.loc[MetricInfo['Metric'] == 'Input_Size','Warn_Cutoff'].values[0]):
            InputSize = ReadDepthHistPlotter(SampleName,_user_df,_bgd_df,1,_figinfo,fig)
            fig = InputSize.Figure

        # Plotting figure 2: Percentage of Reads after Trimming
        if pd.notna(MetricInfo.loc[MetricInfo['Metric'] == 'Percent_PostTrim','Warn_Cutoff'].values[0]):
            TrimmingPercent = TrimmingPlotter(SampleName, _user_df, _bgd_df,2,_figinfo,fig)
            fig = TrimmingPercent.Figure

        # Plotting figure 3: Percentage of Uniquely Aligned Reads
        if pd.notna(MetricInfo.loc[MetricInfo['Metric'] == 'Percent_Uniquely_Aligned','Warn_Cutoff'].values[0]):
            Alignment = AlignmentPlotter(SampleName, _user_df, _bgd_df,3,_figinfo,fig)
            fig = Alignment.Figure

        # Plotting figure 4: Percentage of Reads Mapped to Exons
        if pd.notna(MetricInfo.loc[MetricInfo['Metric'] == 'Percent_Exonic','Warn_Cutoff'].values[0]):
            ExonMapping = ExonMappingPlotter(SampleName, _user_df, _bgd_df,4,_figinfo,fig)
            fig = ExonMapping.Figure

        # Plotting figure 5: Scatter Plot of Number of Ribosomal RNA reads per Uniquely Aligned Reads
        if pd.notna(MetricInfo.loc[MetricInfo['Metric'] == 'Perc_Aligned_Reads_Overlapping_rRNA','Warn_Cutoff'].values[0]):
            fig = plotScatter_rRNA(SampleName, _user_df, _bgd_df, 5,_figinfo,fig)

        # Plotting figure 6: Violin Plot for Contamination - % Adapter Content and % Overrepresented Sequences
        if all(col in _bgd_df.columns for col in ["Percent_Overrepresented_Seq_Trimmed",
                                                  'Percent_Adapter_Content_Trimmed',
                                                  "Percent_Overrepresented_Seq_Untrimmed",
                                                  'Percent_Adapter_Content_Untrimmed']):

            fig = plotViolin_dualAxis(SampleName, _user_df, _bgd_df, 6,_figinfo,fig)

        # Plotting figure 7: Expression Distribution Plot
        if pd.notna(MetricInfo.loc[MetricInfo['Metric'] == 'NumGenes','Warn_Cutoff'].values[0]):
            fig = plotNegBin(SampleName,_user_df,_negBin_df,7,_figinfo,fig)

        # Plotting figure 8: Gene Body Coverage Plot
        if pd.notna(MetricInfo.loc[MetricInfo['Metric'] == 'GBC_KSstats','Warn_Cutoff'].values[0]):
            fig = plotGC(SampleName,_user_df, _gc_df, 8,_figinfo,fig)

        shared_text_kwargs = {'y'                 : .99,
                              'fontsize'          : 6,
                              'verticalalignment' : 'top',
                              'fontweight'        : 'book',
                              'style'             : 'italic'}
        # Add sample info at the top-left corner of the page
        fig.text(s                   = 'Sample : ' + SampleName,
                 x                   = 0.01,
                 horizontalalignment = 'left',
                 **shared_text_kwargs)

        fig.text(s                   = "Batch : " + _user_df.loc[_user_df['Sample'] == SampleName,'Batch'].iloc[0],
                 x                   = 0.99,
                 horizontalalignment = 'right',
                 **shared_text_kwargs)

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

    parser.add_argument("-ref",
                        "--reference-data",
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
    bgd_filename    = args.reference_data if args.reference_data else qry_filename
    cutoff_filename = args.cutoff_data
    fail_alpha      = float(args.failalpha)
    warn_alpha      = float(args.warnalpha)

    print(f"Query Data : {qry_filename}")
    print(f"Output File : {op_folder}")
    print(f"Reference Data : {bgd_filename}")
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
