import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from abc import ABC, abstractmethod
import os
import argparse
matplotlib.use('PDF')
import seaborn
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
#import helper_retroFunction.py
from sklearn.linear_model import LinearRegression
from scipy.stats import norm
from scipy import stats
import matplotlib.offsetbox
from datetime import datetime
from statsmodels.stats.weightstats import ztest
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
from abc import ABC, abstractmethod

# Functions for formatting values....Turn into object?
def fmt_scatter_million(_x,_pos= None):
    return f'{_x:.2E}'

def fmt_million(_x,_pos = None):
    return '{0:.0f}'.format(_x)

def fmt_percent(x,pos = None):
    return f'{x:.2f}%'

def fmt_percent1(x,pos = None):
    return f'{x:.0f}%'

def fmt(_x, _pos):
    return '{0:.1f}'.format(_x)

class UI_adapter:
    """
    The parent class designed to handle the inputs provided by the user, shape it into the correct format, and return 
    meaningful warnings and errors to avoid downstream confusion
    """

    def __init__(self,input_df):
        self.input_df = input_df

    human_readable_names = [] # implemented by subclass
    input_program_names  = [] # implemented by subclass

    def make_mapping_names(self):
            self.mapping_names = dict(zip(self.human_readable_names,self.input_program_names))

    def _validate_column_names(self):
        """Check if the column names of the input DataFrame are included in the input_human_readable_names list"""
        input_columns = set(self.input_df.columns)
        valid_columns = set(self.human_readable_names)
        unrecognized_columns = input_columns.difference(valid_columns)

        if unrecognized_columns:
            raise ValueError(f"Unrecognized columns: {', '.join(unrecognized_columns)}")

    def _fill_missing_values(self):
        """ Fill empty columns if they are missing."""
        self.input_df.fillna(0, inplace = True)

    def _reorder_columns(self):
        """ Reorder columns to fit the format expected by the rest of the program"""
        missing_columns = [col for col in self.input_program_names if col not in self.input_df.columns]
        if missing_columns:
            raise KeyError(f"The following columns are missing from the DataFrame: {missing_columns}")
        self.input_df = self.input_df[self.input_program_names]

    def _change_column_names(self):
        """ Change the user friendly column names to those expected by the internal program"""
        self.input_df.rename(columns=self.mapping_names, inplace=True)

    def _transform_values(self):
        """ Calculate columns that are transformations of the supplied inputs """

        self.input_df["Percent_PostTrim"] =  self.input_df["Num_PostTrim"] / self.input_df["Input_Size"] * 100
        self.input_df["Percent_Uniquely_Aligned"] = self.input_df["Num_Uniquely_Aligned"] / self.input_df["Num_PostTrim"] * 100
        self.input_df["Percent_Exonic"]  = self.input_df["Num_Exonic"] / self.input_df["Num_Uniquely_Aligned"] * 100
        self.input_df['% Aligned Reads Overlapping rRNA'] = self.input_df['Num_Uniquely_Aligned_rRNA'] / self.input_df['Num_Uniquely_Aligned']

    def adapt_input(self):
        # specific methods required by each subclass to be implemented by subclass
        pass

class input_adapter(UI_adapter):
    """
    This object is designed to handle the input provided by the user, shape it into the correct format, and return
    meaningful warnings and errors to avoid downstream confusion

    Input: The loaded pd of the input file
    Output:The standardized file expected by the rest of the program, or an error
    """

    human_readable_names = [
        "Sample",
        "# Sequenced Reads",
        "# Post-trim Reads",
        "Percent_PostTrim", # Calculated by adapter
        "# Uniquely Aligned Reads",
        "% Uniquely Aligned Reads",
        "# Reads Mapped to Exons",
        "Percent_Exonic", # Calculated by adapter
        "# Aligned Reads Overlapping rRNA",
        "% Aligned Reads Overlapping rRNA", # Calculated by adapter
        "% Overrepresented Sequences (Pre-trim)",
        "% Adapter Content (Pre-trim)",
        "% Overrepresented Sequences (Post-trim)",
        "% Adapter Content (Post-trim)",
        "Batch"]
    input_program_names = [
        "Sample",
        "Input_Size",
        "Num_PostTrim",
        "Percent_PostTrim", # Calculated by adapter
        "Num_Uniquely_Aligned",
        "Percent_Uniquely_Aligned", # Calculated by adapter
        "Num_Exonic",
        "Percent_Exonic", # Calculated by adapter
        "Num_Uniquely_Aligned_rRNA",
        "% Aligned Reads Overlapping rRNA", # Calculated by adapter
        "Percent_Overrepresented_Seq_Untrimmed",
        "Percent_Adapter_Content_Untrimmed",
        "Percent_Overrepresented_Seq_Trimmed",
        "Percent_Adapter_Content_Trimmed",
        "Batch"]

    def adapt_input(self):
        self.make_mapping_names()
        self._validate_column_names()
        self._fill_missing_values()
        self._change_column_names()
        self._transform_values()
        self._reorder_columns()

        return self

class manual_cutoff_adapter(UI_adapter):

    """
    This object is designed to handle the manual cutoff file provided by the user, shape it into the correct format, and return
    meaningful warnings and errors to avoid downstream confusion

    Input: The loaded pd of the manual cutoff adapter
    Output:The standardized file expected by the rest of the program, or an error
    """

    human_readable_names = [
        "Cutoff",
        "Sequencing Depth",
        "% of Reads After Trimming",
        "% Uniquely Aligned Reads / Trimmed Reads",
        "% Mapped Reads / Aligned Reads",
        "rRNA Reads / Aligned Reads",
        "% Overrep Sequences (Post-Trim)",
        "% Adapter Content (Post-Trim)",
        "# Detected Genes",
        "Gene Body Coverage (Pval)"]

    input_program_names = [
        "cutoff",
        "Input_Size",
        "Percent_PostTrim",
        "Percent_Uniquely_Aligned",
        "Percent_Exonic",
        "% Aligned Reads Overlapping rRNA",
        "Percent_Overrepresented_Seq_Trimmed",
        "Percent_Adapter_Content_Trimmed",
        "NumGenes",
        "GC_cutoff"]

    def _transpose_df(self):
        """transposes the df with additional operations to get the right columns names"""
        self.input_df = self.input_df.transpose()
        self.input_df.columns = self.input_df.iloc[0]
        self.input_df = self.input_df.drop(self.input_df.index[0])

    def adapt_input(self):
        self.make_mapping_names()
        self._validate_column_names()
        self._change_column_names()
        self._reorder_columns()
        self._transpose_df()

def add_warn_fail_markers(_figinfo,ax,cutoff_key):

    def add_vert_marker(ax,cutoff,plot_scalar,clr,txt):

        ax.plot(cutoff,
                ax.get_ylim()[1] - (18 * plot_scalar),
                marker = 'v',
                ms     = 0.8,
                c      = clr)
        ax.text(cutoff,
                ax.get_ylim()[1] - (13 * plot_scalar),
                txt,
                fontsize            = 4,
                color               = clr,
                horizontalalignment ='center')
        return ax

    plot_scalar =  get_axis_range(ax.get_ylim()) / 100

    ax = add_vert_marker(ax          = ax,
                         cutoff      = _figinfo["_fail_cutoffs"][cutoff_key],
                         plot_scalar = plot_scalar,
                         clr         = _figinfo["_fail_color"],
                         txt         = "Fail")
    ax = add_vert_marker(ax          = ax,
                         cutoff      = _figinfo["_warn_cutoffs"][cutoff_key],
                         plot_scalar = plot_scalar,
                         clr         = _figinfo["_warn_color"],
                         txt         = "Warn")

    return ax

# This function calculates whether the current sample label orientation needs 
# to be adjusted and returns the required vars
def adjust_flag(_ax,_current_sample,_lib_mean,Formatter= fmt_scatter_million):

    labelpadding = (get_axis_range(_ax.get_xlim())/ 50)

    curr_sample_kwargs = {'y'        : (_ax.get_ylim()[1] / 3.5),
                          's'        : Formatter(_current_sample),
                          'fontsize' : 3,
                          'zorder'   : 2}
    lib_kwargs = curr_sample_kwargs.copy()
    lib_kwargs.update({'s' : Formatter(_lib_mean)})
    if _current_sample >=  _lib_mean:
        # plot the current sample line
        _ax.text(x        = _current_sample + (labelpadding / 2),
                 rotation = 270,
                 **curr_sample_kwargs)
        # plot the library mean line
        _ax.text(x        = _lib_mean - labelpadding,# this needs to be halved to look good
                 rotation = 90,
                 **lib_kwargs)
    else:
        # plot the current sample line
        _ax.text(x        = _current_sample - labelpadding,
                 rotation = 90,
                 **curr_sample_kwargs)
        # plot the library mean line
        _ax.text(x        = _lib_mean + (labelpadding / 2),
                 rotation = 270,
                 **lib_kwargs)
    return _ax

# set the legend for figures with warn / fail IE 1_6)
def legend_setup_1_6(_ax,_line1,_line2,_figinfo,cutoff_key,_loc,formatter = fmt_million):
    _fail_label = mpatches.Patch(color = _figinfo["_fail_color"],
                                 label = 'Fail Cutoff')
    _warn_label=  mpatches.Patch(color = _figinfo["_warn_color"],
                                 label = 'Warn Cutoff')
    _ax.legend(handles   = [_line1,
                            _line2,
                            _fail_label,
                            _warn_label],
               labels    = ["Current Sample",
                            "Batch Mean",
                            "Fail (" + formatter(_figinfo["_fail_cutoffs"][cutoff_key]) + ")",
                            "Warn (" + formatter(_figinfo["_warn_cutoffs"][cutoff_key]) + ")"],
                loc      = _loc,
                frameon  = False,
                fontsize = _figinfo["_legend_size"])
    return _ax

# the goal of this function is to set which axes of the regular axis and the Kernel density axis 
def mk_axes(_plt_ax,_kd_ax = None):

    for label in _plt_ax.get_xticklabels() + _plt_ax.get_yticklabels():
        label.set_fontsize(4)

    # Configure plot axis appearance
    for spine in ['top', 'right']:
        _plt_ax.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']:
        _plt_ax.spines[spine].set_visible(True)
        _plt_ax.spines[spine].set_color('black')
        _plt_ax.spines[spine].set_linewidth(0.55)
    _plt_ax.set_facecolor('white')

    if _kd_ax is not None:
        # Make kernel density axis invisible
        _kd_ax.yaxis.set_visible(False)
        _kd_ax.xaxis.set_visible(False)
        for spine in ['top', 'right', 'bottom', 'left']:
            _kd_ax.spines[spine].set_visible(False)
        return _plt_ax, _kd_ax

    return _plt_ax

# the goal of this function is to determine if a plot needs a fail or warn box and then call the according
# helper function 
def needs_fail_or_warn(ax,current_sample,_figinfo,cutoff_key,higher_lower):

    class FlagInserter:
        def __init__(self, ax=None):
            self.ax = ax or plt.gca()
            self.cutoff_warn = _figinfo["_warn_cutoffs"][cutoff_key]
            self.cutoff_fail = _figinfo["_fail_cutoffs"][cutoff_key]

        def make_flag(self, flag_type):
            if flag_type == "fail":
                text = "FAILURE"
                background_color = 'tomato'
                text_color = 'yellow'
            elif flag_type == "warn":
                text = "WARNING"
                background_color = 'yellow'
                text_color = 'red'
            else:
                raise ValueError("Invalid flag type")

            anch_text = matplotlib.offsetbox.AnchoredText(text,
                                                          pad            = 0.0001,
                                                          borderpad      = 1,
                                                          loc            = 3,
                                                          prop           = dict(size            = 4,
                                                                                snap            = True,
                                                                                backgroundcolor = background_color,
                                                                                color           = text_color,
                                                                                fontweight      = 'roman',
                                                                                fontfamily      = 'serif',
                                                                                fontstretch     = 'expanded'),
                                                          frameon        = True,
                                                          bbox_to_anchor = (0, 1),
                                                          bbox_transform = self.ax.transAxes)

            self.ax.add_artist(anch_text)
            return None

    def insert_flag(ax, cutoff, flag_func):
        if higher_lower == "lower" and current_sample <= cutoff  or higher_lower == "upper" and current_sample >= cutoff:
            flag_func(ax)

    flag_inserter = FlagInserter()

    if higher_lower == "lower":
        if current_sample <= flag_inserter.cutoff_fail:
            insert_flag(ax, flag_inserter.cutoff_fail, lambda ax: flag_inserter.make_flag("fail"))
        elif current_sample <= flag_inserter.cutoff_warn:
            insert_flag(ax, flag_inserter.cutoff_warn, lambda ax: flag_inserter.make_flag("warn"))
    elif higher_lower == "upper":
        if current_sample >= flag_inserter.cutoff_fail:
            insert_flag(ax, flag_inserter.cutoff_fail, lambda ax: flag_inserter.make_flag("fail"))
        elif current_sample >= flag_inserter.cutoff_warn:
            insert_flag(ax, flag_inserter.cutoff_warn, lambda ax: flag_inserter.make_flag("warn"))

    return ax

# The goal of this function is to return the upper or/and lower bound of a ci given a vec
def get_ci_bound(vec, alpha, upper_lower="both"):
    mean = np.mean(vec)
    scale = stats.tstd(vec)
    lower, upper = stats.norm.interval(alpha= 1 - alpha, loc=mean, scale=scale)

    if upper_lower == "upper":
        return upper
    elif upper_lower == "lower":
        return lower
    elif upper_lower == "both":
        return lower, upper
    else:
        raise ValueError("Invalid value for 'upper_lower'. Must be 'upper', 'lower', or 'both'.")

# Simple bootstrap calculation
def CalcBootstrapBound(vec,upper_lower,alpha,num_resamples = 1000):

    """
    n = len(vec)
    bootstrap_medians = []
    for _ in range(num_resamples):
        resampled_data = np.random.choice(vec,
                                          size = n,
                                          replace = True)
        median = np.mean(resampled_data)
        bootstrap_medians.append(median)

    bootstrapped_values = bootstrap_medians

    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100

    lower_bound = np.percentile(bootstrapped_values, lower_percentile)
    upper_bound = np.percentile(bootstrapped_values, upper_percentile)

    if upper_lower == "upper":
        return upper_bound
    if upper_lower == "lower":
        return lower_bound
    """
    vec = np.asarray(vec)
    bootstrap_mean = bs.bootstrap(vec,
                                  stat_func = bs_stats.median).value
    bootstrap_std  = bs.bootstrap(vec,
                                  stat_func = bs_stats.std).value
    conf_size = norm.ppf(alpha) * bootstrap_std
    if upper_lower == "upper":
        cutoff = bootstrap_mean - conf_size
    if upper_lower == "lower":
        cutoff = bootstrap_mean + conf_size

    return np.float64(cutoff)

class CutoffCalculator:
    def __init__(self, bgd_df, alph):
        self.bgd_df = bgd_df
        self.onesided_alph = 2*alph
        self.cutoffs_dict = {}
        self.cutoffs_dict["_alpha"] = alph

    def __call__(self):

        for metric in ["Input_Size","Percent_PostTrim","Percent_Uniquely_Aligned","Percent_Exonic"]:
            self.calculate_cutoff(metric,'lower')

        for metric in ["% Aligned Reads Overlapping rRNA","Percent_Overrepresented_Seq_Trimmed","Percent_Adapter_Content_Trimmed"]:
            self.calculate_cutoff(metric, "upper")

        return self.cutoffs_dict

    def calculate_cutoff(self,column,upper_lower):
        vec = np.array(self.bgd_df.loc[:, column])
        bootstrap_mean = bs.bootstrap(vec,
                                      stat_func = bs_stats.mean).value
        bootstrap_std  = bs.bootstrap(vec,
                                      stat_func = bs_stats.std).value
        conf_size = norm.ppf(self.onesided_alph) * bootstrap_std
        if upper_lower == "upper":
            self.cutoffs_dict[column] = bootstrap_mean - conf_size
        if upper_lower == "lower":
            self.cutoffs_dict[column] = bootstrap_mean + conf_size


def gen_cutoffs(bgd_df, alph):
    calculator = CutoffCalculator(bgd_df, alph)
    return calculator()

def set_ticks(_ax,_tick_size):

    shared_tick_kwarggs = {'which'     : 'both',
                          'length'    : 1,
                          'width'     : .5,
                          'labelsize' : _tick_size,
                          'direction' : 'out',
                          'pad'       : 2}
    _ax.tick_params(axis        = 'x',
                    labelbottom = True,
                    bottom      = True,
                    **shared_tick_kwarggs)
    _ax.tick_params(axis      = 'y',
                    labelleft = True,
                    left      = True,
                    **shared_tick_kwarggs)

    return _ax

# takes two vectors, finds the maximum and minimum, and creates the number of bins specified
def make_bins(vec1,vec2,num_bins):

    plt_min,plt_max = pd.concat([vec1,vec2]).agg(['min','max'])
    bins = np.linspace(plt_min,
                       plt_max,
                       num_bins + 1)

    return bins

def mkTitlePage(_figinfo):

    """
    Make QC heatmap data
    The Goal of this function is to calculate a matrix that can be used with matplotlibs heatmap
    function which contains whether a sample passed, failed, or was warned for each test
    """


    def mk_cutoff_descript(descriptor,cutoff_set):
        print('the cutoff set is',cutoff_set)

        GC_cutoff_text = (f"Gene Body Coverage K-S Stat : >{cutoff_set['GC_cutoff']:.3f} | "
                          if cutoff_set['GC_cutoff'] else "")
        Hist_cutoff_text =  (f"# Detected Genes : < {int(cutoff_set['NumGenes'])}"
                             if cutoff_set['NumGenes'] else "")
        cutoff_descriptor =(f"{descriptor}: Default alpha : {cutoff_set['_alpha']} | "
                            f"Sequencing Depth : < {int(cutoff_set['Input_Size'])} | "
                            f"Trimming : < {cutoff_set['Percent_PostTrim']:.3f}% | "
                            f"Alignment : < {cutoff_set['Percent_Uniquely_Aligned']:.3f}% | "
                            f"Gene Exon Mapping : < {cutoff_set['Percent_Exonic']:.3f}% | "
                            f"Ribosomal RNA : > {cutoff_set['% Aligned Reads Overlapping rRNA']:.3f}% |\n "
                            f"Adapter Contamination : >{cutoff_set['Percent_Adapter_Content_Trimmed']:.3f}% | "
                            f"Overrep. Seq Contamination : >{cutoff_set['Percent_Overrepresented_Seq_Trimmed']:.3f}% | "
                            f"{GC_cutoff_text}"
                            f"{Hist_cutoff_text}")
        return cutoff_descriptor

    fig = plt.figure()

    # add text
    fig.text(.5,
             .965,
             "QC Plotter Input Summary",
             ha       = 'center',
             va       = 'top',
             fontsize = 14)
    fig.text(.005,
             .985,
             datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
             fontsize = 4)
    fig.text(.005,
             .97,
             "Version 2.0",
             fontsize = 4)
    fig.text(.5,
             .5,
             (f"Input table = {_figinfo['_ip_filename']} \n"
              f"Output location =    {_figinfo['_op_filename']}\n"
              f"Background table =    {_figinfo['_bgd_filename']}\n"
              f"Gene Body Coverage file =    {_figinfo['_gc_file']}\n"
              f"Raw count table  = {_figinfo['_hist_file']}\n"
              f"Cutoff file =     {_figinfo['_cutoff_filename']}"),
             fontsize = 6,
             ha       = "center",
             va       = "center")

    # show cutoffs
    warn_descript = mk_cutoff_descript("Warn cutoffs",_figinfo["_warn_cutoffs"]) 
    fail_descript = mk_cutoff_descript("Fail cutoffs",_figinfo["_fail_cutoffs"]) 


    shared_text_characteristics = {'fontsize'   : 5,
                                   'ha'         : 'center',
                                   'va'         : 'top',
                                   'fontweight' : 'book',
                                   'style'      : 'italic'}
    fig.text(s          = warn_descript,
             x          = .5,
             y          = .1,
             **shared_text_characteristics)

    fig.text(s          = fail_descript,
             x          = .5,
             y          = .05,
             **shared_text_characteristics)

    return fig

# Used for calculating 
class MetricStatusStrategy(ABC):
    @abstractmethod
    def compute_status(self, test_value, warn_value, fail_value):
        pass

class LowerStatusStrategy(MetricStatusStrategy):
    def compute_status(self, test_value, warn_value, fail_value):
        if test_value <= warn_value:
            return 1 if test_value <= fail_value else 0.5
        return 0

class UpperStatusStrategy(MetricStatusStrategy):
    def compute_status(self, test_value, warn_value, fail_value):
        return LowerStatusStrategy().compute_status(-test_value, -warn_value, -fail_value)

def mkQC_heatmap_data(_userDf, _figinfo):
    lower_status_strategy = LowerStatusStrategy()
    upper_status_strategy = UpperStatusStrategy()

    strategies = [
        ("Input_Size", "Input_Size", lower_status_strategy),
        ("Percent_PostTrim", "Percent_PostTrim", lower_status_strategy),
        ("Percent_Uniquely_Aligned","Percent_Uniquely_Aligned",lower_status_strategy),
        ("Percent_Exonic","Percent_Exonic",lower_status_strategy),
        ("% Aligned Reads Overlapping rRNA","% Aligned Reads Overlapping rRNA",upper_status_strategy),
        ("Percent_Overrepresented_Seq_Trimmed","Percent_Overrepresented_Seq_Trimmed",upper_status_strategy),
        ("Percent_Adapter_Content_Trimmed","Percent_Adapter_Content_Trimmed",upper_status_strategy)]
    if _figinfo["_fail_cutoffs"]['NumGenes'] is not None:
        strategies.append(("NumGenes","NumGenes",lower_status_strategy))

    if _figinfo["_fail_cutoffs"]['GC_cutoff'] is not None:
        strategies.append(("GBC_KSstats","GC_cutoff",upper_status_strategy))

    _htmat = np.zeros((len(_userDf), 9))

    for _tuple in _userDf.itertuples():
        for i, (column, key, strategy) in enumerate(strategies):

            test_value = _userDf.iloc[_tuple.Index][column]
            warn_value = _figinfo["_warn_cutoffs"][key]
            fail_value = _figinfo["_fail_cutoffs"][key]

            _htmat[_tuple.Index, i] = strategy.compute_status(test_value, warn_value, fail_value)

    return _htmat

def mkQC_heatmap(heatmap_data):

    """
    mkQC_heatmap: This function generates the summary heatmap
    take input from mkQC_heatmap_data as input
    """

    # get how many rows there are
    numrows, numcols = heatmap_data.shape
    cellwidth = .25
    cellheight = .1

    fig_width = numcols * cellwidth
    fig_height= numrows * cellheight

    page_height = 5
    page_width = 6.3
    colors = ["grey","goldenrod","red"]
    cm = matplotlib.colors.ListedColormap(colors)
    sample_names = heatmap_data.Sample
    heatmap_data = heatmap_data.drop("Sample",axis = 1)
    fig2,ax = plt.subplots(figsize=(page_width, page_height))
    fig2.text(s        = "Summary of QC Metrics",
              x        = .5,
              y        = .9,
              fontsize = 10,
              ha       = 'center')
    seaborn.heatmap(heatmap_data.values,
                    ax          = ax,
                    xticklabels = ["Sequencing Depth","Trimming","Alignment","Exon Mapping","Ribosomal RNA",
                                   "Sequence Contamination (Overrep)","Sequence Contamination (Adapter)",
                                   "# Detected Genes","Gene Body Coverage"],
                    yticklabels = sample_names,
                    cmap        = cm,
                    cbar_kws    = {'aspect' : 2,
                                   'shrink' : 2,
                                   'pad'    : .025})
    # change y-axis tick label font size

    width_padding = (1 - fig_width / page_width) / 2
    height_padding = (1 - fig_height / page_height) / 2

    fig2.subplots_adjust(left   =  width_padding,
                         right  = 1 -  width_padding,
                         top    = 1 - height_padding,
                         bottom = height_padding)

    ax.set_yticklabels(ax.get_yticklabels(), fontsize = 5)
    # change x-axis tick label font size
    ax.set_xticklabels(ax.get_xticklabels(), fontsize = 5)
    plt.yticks(rotation = 0)
    plt.xticks(rotation = 90)
    # Get the Colorbar object from the heatmap
    cbar = ax.collections[0].colorbar
    cbar.ax.set_aspect(.5)
    cbar.ax.set_anchor("N")
    # Change the ticks on the colorbar
    cbar.set_ticks([0, 0.5, 1])
    ax.collections[0].colorbar.ax.tick_params(labelsize=5)
    cbar.set_ticklabels(['Passed', 'Warned', 'Failed'])

#    plt.subplots_adjust(left=0.3, bottom=0.3, right=0.7, top=0.8)

    return fig2


class AbstractHistPlotter(ABC):

    # Default values for subclass-specific attributes. Overwritten by concrete subclasses
    VarName   = None
    CutoffKey = None
    Formatter = None
    PlotTitle = None
    XAxisTitle= None

    def __init__(self,SampleName, _user_df, _background_df, _position,_figinfo,_figure=None):
        self.SampleName = SampleName #Which sample"
        self.Position= _position#"Where to put the plot"
        self.BgdDf   = _background_df
        self.UserDf  = _user_df
        self.FigInfo = _figinfo
        self.Figure  = _figure  #"The Sample Sheet being added to"

    def InitDependentFields(self):
        self.UserVals = self.UserDf[self.VarName]
        self.BgdVals  = self.BgdDf[self.VarName] #"the background to compare the sample to"

    def AddHist(self):
        _bins = make_bins(self.BgdVals,
                          self.UserVals,
                          self.FigInfo["_bin_num"])

        _current_sample = self.UserDf[self.UserDf['Sample'] == self.SampleName].iloc[0]
        _current_value = _current_sample[self.VarName]
        _lib_mean = self.UserDf.loc[self.UserDf['Batch'] == _current_sample['Batch'],self.VarName].mean()
        warn_limit = self.FigInfo['_warn_cutoffs'][self.CutoffKey]
        fail_limit = self.FigInfo['_fail_cutoffs'][self.CutoffKey]

        axis = self.Figure.add_subplot(self.FigInfo["_subplot_rows"],
                                       2,
                                       self.Position)

        sns.histplot(self.BgdVals,
                     bins      = _bins,
                     ax        = axis,
                     color     = 'lightgray',
                     edgecolor = "lightgray")

        axis1 = axis.twinx()
        sns.kdeplot(self.BgdVals,
                    ax        = axis1,
                    color     = 'black',
                    lw        = 0.5,
                    bw_adjust = .5)

        # set limits
        _xmin = min(self.BgdVals.min(),
                    _lib_mean,
                    warn_limit,
                    fail_limit,
                    _current_value)
        _xmax = max(self.BgdVals.max(),
                    _lib_mean,
                    warn_limit,
                    fail_limit,
                    _current_value)
        cushion = .05 * (_xmax - _xmin)
        axis.set_xlim(_xmin - cushion,_xmax + cushion)
        axis = set_ticks(axis, self.FigInfo["_tick_size"])

        axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
        axis.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
        axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(self.Formatter))

        axis.set_title(self.PlotTitle,
                       fontsize = self.FigInfo["_title_size"])

        axis.set_xlabel(self.XAxisTitle,
                        labelpad = 1,
                        fontsize = self.FigInfo["_label_size"])

        axis.set_ylabel('# Samples',
                        labelpad = 2,
                        fontsize = self.FigInfo["_label_size"])

        axis = adjust_flag(axis,
                           _current_value,
                           _lib_mean,
                           self.Formatter)

        ### Adding cutoff markers
        axis = add_warn_fail_markers(self.FigInfo,
                                     axis,
                                     self.CutoffKey)

        # Current Sample Line and Label
        SampleLine = axis.axvline(x         = _current_value,
                                  color     = self.FigInfo["_curr_sample_color"],
                                  linestyle = '-',
                                  linewidth = 0.5,
                                  label     = self.Formatter(_current_value))

        # Current Library Mean Line and Label
        BgdLine = axis.axvline(x        = _lib_mean,
                               color     = 'indigo',
                               linestyle = '--',
                               linewidth = 0.5,
                               label     = self.Formatter(_lib_mean))

        # set up axes
        axis = legend_setup_1_6(_ax        = axis,
                                _line1     = SampleLine,
                                _line2     = BgdLine,
                                _figinfo   = self.FigInfo,
                                _loc       = 'upper left',
                                cutoff_key = self.CutoffKey,
                                formatter  = self.Formatter)

        #set axes to be visible or not
        axis,axis1 =  mk_axes(axis,axis1)

        axis = needs_fail_or_warn(ax             = axis,
                                  current_sample = _current_value,
                                  _figinfo       = self.FigInfo,
                                  cutoff_key     = self.CutoffKey,
                                  higher_lower   = "lower")

class ReadDepthHistPlotter(AbstractHistPlotter):

    def __init__(self, _ip_tuple, _user_df, _background_df, _position, _figinfo, _figure=None):
        self.VarName    =  "Input_Size"
        self.CutoffKey  = "Input_Size"
        self.Formatter  = fmt_scatter_million
        self.PlotTitle  = "Sequencing Depth"
        self.XAxisTitle = "# Sequenced Reads"
        super().__init__(_ip_tuple, _user_df, _background_df, _position, _figinfo, _figure)
        self.InitDependentFields()
        self.AddHist()

class TrimmingPlotter(AbstractHistPlotter):

    def __init__(self, _ip_tuple, _user_df, _background_df, _position, _figinfo, _figure=None):
        self.VarName   = "Percent_PostTrim"
        self.CutoffKey = "Percent_PostTrim"
        self.Formatter = fmt_percent
        self.PlotTitle = "Trimming"
        self.XAxisTitle = "% Post-trim Reads"
        super().__init__(_ip_tuple, _user_df, _background_df, _position, _figinfo, _figure)
        self.InitDependentFields()
        self.AddHist()

class AlignmentPlotter(AbstractHistPlotter):

    def __init__(self, _ip_tuple, _user_df, _background_df, _position, _figinfo, _figure=None):
        self.VarName  = "Percent_Uniquely_Aligned"
        self.CutoffKey= "Percent_Uniquely_Aligned"
        self.Formatter= fmt_percent
        self.PlotTitle= "Alignment"
        self.XAxisTitle= "% Uniquely Aligned Reads"
        super().__init__(_ip_tuple, _user_df, _background_df, _position, _figinfo, _figure)
        self.InitDependentFields()
        self.AddHist()

class ExonMappingPlotter(AbstractHistPlotter):

    def __init__(self, _ip_tuple, _user_df, _background_df, _position, _figinfo, _figure=None):
        self.VarName  = "Percent_Exonic"
        self.CutoffKey= "Percent_Exonic"
        self.Formatter= fmt_percent
        self.PlotTitle= "Exon Mapping"
        self.XAxisTitle= "% Mapped to Exons / Aligned Reads"
        super().__init__(_ip_tuple, _user_df, _background_df, _position, _figinfo, _figure)
        self.InitDependentFields()
        self.AddHist()

#### Plot 5: rRNA Scatter ####
def plotScatter_rRNA(SampleName, _userDf, _background_df, _pos,_figinfo,_f=None):

    _ax = plt.subplot(_figinfo["_subplot_rows"],
                      2,
                      _pos)

    _plotter_df = pd.concat([_background_df,
                             _userDf],
                            sort = True)

    # Assign color for current project's library (all samples in the current project)
    current_batch = _userDf[_userDf['Sample'] == SampleName].iloc[0]['Batch']

    _plotter_df["scatter_color"] = np.where(_plotter_df["Batch"] == current_batch,
                                            "indigo",
                                            "lightgray")

    # Assign separate color for current sample on each page
    _plotter_df.loc[_plotter_df["Sample"] == SampleName,
                    "scatter_color"] = _figinfo["_curr_sample_color"]
    _intupdf = _plotter_df.loc[_plotter_df["Sample"] == SampleName]

    ## Regression line (gradient slope)
    X = _background_df.loc[:, "Num_Uniquely_Aligned"].values.reshape(-1, 1)
    Y = _background_df.loc[:, "Num_Uniquely_Aligned_rRNA"].values.reshape(-1, 1)
    linear_regressor = LinearRegression()
    linear_regressor.fit(X, Y)
    Y_pred = linear_regressor.predict(X)
    _slope_current = linear_regressor.coef_[0][0]

    # plot the regression for all points.
    _ax.scatter(x = _plotter_df['Num_Uniquely_Aligned'],
                y = _plotter_df['Num_Uniquely_Aligned_rRNA'],
                s = 0.8,
                c = _plotter_df["scatter_color"])

    #separate scatter call for the sample so it can have a unique size and shape
    _ax.scatter(x      =_intupdf['Num_Uniquely_Aligned'],
                y      =_intupdf['Num_Uniquely_Aligned_rRNA'],
                marker = "*",
                s      =20,
                c      =_intupdf["scatter_color"])

    _ax.set_title("Ribosomal RNA",
                  fontsize = _figinfo["_title_size"])
    _ax = set_ticks(_ax,
                    _figinfo["_tick_size"])

    _ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_scatter_million))
    _ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_scatter_million))
    _ax.set_xlabel("# Uniquely Aligned Reads",
                   fontsize= _figinfo["_label_size"],
                   labelpad =2)

    _ax.set_ylabel("# Aligned rRNA Reads",
                   fontsize= _figinfo["_label_size"],
                   labelpad= 2)

    xmin, xmax = _ax.get_xlim()
    line_x0 = 0
    line_y0 = 0

    line_y1_warn = _figinfo["_warn_cutoffs"]["% Aligned Reads Overlapping rRNA"] * (xmax - line_x0) + line_y0
    line_y1_fail = _figinfo["_fail_cutoffs"]["% Aligned Reads Overlapping rRNA"] * (xmax - line_x0) + line_y0

    mean_line_kwargs = {'linestyle' : '--',
                        'linewidth' : .7}

    # Plot the regression line
    _ax.plot(X,
             Y_pred,
             c         = 'black',
             **mean_line_kwargs)
    _ax.plot([line_x0,
              xmax],
             [line_y0,
              line_y1_warn],
             c         =_figinfo["_warn_color"],
             label     = "Warn",
             **mean_line_kwargs)
    _ax.plot([line_x0,
              xmax],
             [line_y0,
              line_y1_fail],
             c         = _figinfo["_fail_color"],
             label     = "Fail",
             **mean_line_kwargs)
    # Set axes margins for padding on both axes

    _ax.set_aspect('auto',
                   adjustable = 'box',
                   anchor     = 'SW')

    _curr_lib = matplotlib.lines.Line2D([0],
                                        [0],
                                        color           = 'w',
                                        markerfacecolor = 'indigo',
                                        marker          = 'o',
                                        markersize      = 3.5)
    _curr_samp = matplotlib.lines.Line2D([0],
                                         [0],
                                         color           = 'w',
                                         markerfacecolor = _figinfo["_curr_sample_color"],
                                         marker          = '*',
                                         markersize      = 6)

    _mean_label = mpatches.Patch(color = 'black',
                                 label = 'Mean Slope')
    _fail_label = mpatches.Patch(color = _figinfo["_fail_color"],
                                 label = 'Fail Cutoff')
    _warn_label = mpatches.Patch(color = _figinfo["_warn_color"],
                                 label = 'Warn Cutoff')

    _ax.legend(handles   = [_curr_samp,
                            _curr_lib,
                            _fail_label,
                            _warn_label,
                            _mean_label],
               labels    =  ["Current Sample",
                             "Batch Samples",
                             f"Fail ({_figinfo['_fail_cutoffs']['% Aligned Reads Overlapping rRNA']:.0%})",
                             f"Warn ({_figinfo['_warn_cutoffs']['% Aligned Reads Overlapping rRNA']:.0%})",
                             f"Mean rRNA/Aligned Reads ({_slope_current:.0%})"],
                loc      = 'upper left',
                frameon  = False,
                fontsize = _figinfo["_legend_size"])

    _ax = mk_axes(_ax)
    _ax = needs_fail_or_warn(_ax,
                             _slope_current,
                             _figinfo,
                             "% Aligned Reads Overlapping rRNA",
                             "upper")

    return _f

#### Plot 6: Sequence Contamination - Violin Plot ####

def plotViolin_dualAxis(SampleName, _userDf, _background_df, _position,_figinfo,_f=None):

    # for plotting the individual composite functions
    def SingleViolin(_axis,_current_sample,_background_df,OverrepName,AdaptName):

        _contaminant_df = _background_df[[OverrepName,
                                          AdaptName]]
        _contaminant_df.columns = ["Overrepresented",
                                   "Adapter"]

        _contaminant_melt  = pd.melt(_contaminant_df,
                                     var_name   = "Contamination_Metric",
                                     value_name = "Percent")

        _current_overrep_untrim = _current_sample[OverrepName].iloc[0]
        _current_adapter_untrim = _current_sample[AdaptName].iloc[0]

        # Format Axes
        _axis.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5))
        _axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_percent))
        _axis = mk_axes(_axis)
        _axis = set_ticks(_axis,
                          _figinfo['_tick_size'])

        # Add violin plot
        sns.violinplot(x         = "Percent",
                       y         = "Contamination_Metric",
                       data      = _contaminant_melt,
                       palette   = _contaminant_pal,
                       inner     = None,
                       ax        = _axis,
                       linewidth = 0.3,
                       orient    = "h",
                       scale     = "count")
        _axis.set_xlabel('')
        _axis.set_ylabel('')

        # Add  4 lines one for the actual value and then the mean for all metrics
        # Define how the lines for each item should be visualized 
        _linekwargs_overrep = {'linestyle' : '-',
                               'linewidth' : .35,
                               'ymin'      : .5,
                               'ymax'      : .95}
        _linekwargs_adapter = _linekwargs_overrep.copy()
        _linekwargs_adapter.update({'ymin'      : .05,
                                   'ymax'      : .45})

        current_batch = _current_sample['Batch'].iloc[0]

        _mean_overrep_untrim = _background_df.loc[_background_df['Batch'] == current_batch,  OverrepName].mean()
        _mean_adapter_untrim = _background_df.loc[_background_df['Batch'] == current_batch, AdaptName].mean()

        _line_overrep = _axis.axvline(x        = _current_overrep_untrim,
                                      color    = _figinfo["_curr_sample_color"],
                                      label    = '{:.2f}%'.format(_current_overrep_untrim),
                                      **_linekwargs_overrep)

        _line_mean_overrep = _axis.axvline(x        = _mean_overrep_untrim,
                                           color    = 'indigo',
                                           label    = '{:.2f}%'.format(_mean_overrep_untrim),
                                           **_linekwargs_overrep)

        _line_adapter = _axis.axvline(x        = _current_adapter_untrim,
                                      color    = _figinfo["_curr_sample_color"],
                                      label    = '{:.2f}%'.format(_current_adapter_untrim),
                                      **_linekwargs_adapter)

        _line_mean_adapter = _axis.axvline(x        = _mean_adapter_untrim,
                                           color    = "indigo",
                                           label    = '{:.2f}%'.format(_mean_adapter_untrim),
                                           **_linekwargs_adapter)

        _axis.legend(handles  = [_line_overrep,
                                 _line_mean_overrep],
                     labels   = ["Current Sample", "Batch Mean"],
                     loc      = 'upper right',
                     frameon  = False,
                     ncol     = 1,
                     fontsize = _figinfo["_legend_size"])
        # Set appealing limits
        _axis.margins(x = .05)

        return None # No need to return anything because the _axis is modified inline

    # Define color palette
    _contaminant_pal = {"Overrepresented": "lightgray",
                        "Adapter"        : "gray"}

    # Specify the locations of the individual stuff
    _gridsp = matplotlib.gridspec.GridSpec(_figinfo["_subplot_rows"]*2,
                                           2,
                                           figure=_f)

    _current_sample = _userDf.loc[_userDf['Sample'] == SampleName]

    _axis  = _f.add_subplot(_gridsp[4, 1:])

    # Create the composing violin plots
    SingleViolin(_axis           = _axis,
                 _current_sample = _current_sample,
                 _background_df  = _background_df,
                 OverrepName     = 'Percent_Overrepresented_Seq_Untrimmed',
                 AdaptName       = 'Percent_Adapter_Content_Untrimmed')

    # Add warn and fail flags
    needs_fail_or_warn(ax             = _axis,
                       current_sample = _current_sample['Percent_Overrepresented_Seq_Trimmed'].iloc[0],
                       _figinfo       = _figinfo,
                       cutoff_key     = "Percent_Overrepresented_Seq_Trimmed",
                       higher_lower   = "upper")

    needs_fail_or_warn(ax             = _axis,
                       current_sample = _current_sample['Percent_Adapter_Content_Trimmed'].iloc[0],
                       _figinfo       = _figinfo,
                       cutoff_key     = "Percent_Adapter_Content_Trimmed",
                       higher_lower   = "upper")

    _axis2 = _f.add_subplot(_gridsp[5, 1:])
    SingleViolin(_axis           = _axis2,
                 _current_sample = _current_sample,
                 _background_df  = _background_df,
                 OverrepName     = 'Percent_Overrepresented_Seq_Trimmed',
                 AdaptName       = 'Percent_Adapter_Content_Trimmed')

    # Format axes and titles
    _axis2.set_xlabel(xlabel   = "% of Reads",
                      fontsize = _figinfo["_label_size"],
                      labelpad = 0.5)

    # stuff at the end`
    _axis.set_title(label    = "Sequence Contamination",
                    fontsize = _figinfo["_title_size"],
                    pad      = 0)

    _axis.set_yticklabels(['Untrimmed \nOverrepresented',
                           'Untrimmed \nAdapter'])
    _axis2.set_yticklabels(['Trimmed \nOverrepresented',
                            'Trimmed \nAdapter'])

    _axis2.legend().remove()

    ### Adding cutoff markers
    _axis3 = _axis2.twinx()
    _axis2,_axis3 = mk_axes(_axis2,_axis3)
    _axis3.yaxis.set_ticks([])
    _axis3.xaxis.label.set_visible(False)
    _axis3.set_ylim(_axis2.get_ylim()[0], _axis2.get_ylim()[1])

    _markers = [ # overrepresented
    (_axis3, _figinfo["_fail_cutoffs"]["Percent_Overrepresented_Seq_Trimmed"], -0.4, 'Fail', _figinfo["_fail_color"]),
    (_axis3, _figinfo["_warn_cutoffs"]["Percent_Overrepresented_Seq_Trimmed"], -0.6,'Warn', _figinfo["_warn_color"]),
    # adapter 
    (_axis3, _figinfo["_fail_cutoffs"]["Percent_Adapter_Content_Trimmed"], .75, 'Fail', _figinfo["_fail_color"]),
    (_axis3, _figinfo["_warn_cutoffs"]["Percent_Adapter_Content_Trimmed"], .85, 'Warn', _figinfo["_warn_color"])]

    for _axs, _cutoff, yloc,label, color in _markers:
        _axs.plot(_cutoff,
                  yloc,
                  marker  = 'v',
                  ms      = 1,
                  c       = color,
                  clip_on = False)
        _axs.text(x        = _cutoff,
                  y        = yloc - .1 ,
                  s        = label,
                  fontsize = _figinfo["_tick_size"],
                  color    = color,
                  ha       = 'center')

    plt.subplots_adjust(hspace = 0)

    return _f

def EstimateDeviances(data_df):
    mean_sample = data_df.mean(axis = 1) # compute mean across columns for each row

    deviances = data_df.sub(mean_sample,
                            axis = 0).abs()
    return deviances

def calculate_distribution_diff_pvals(data_df, n_bootstraps = 1000):

    deviances = EstimateDeviances(data_df)
    deviances = data_df.sub(mean_sample,
                            axis = 0).abs()

    def bootstrap_deviances(deviances, n_bootstraps = 1000):

        bootstrap_distributions = []

        for _ in range(n_bootstraps):

            bootstrap_sample = deviances.sample(frac    = 1,
                                                replace = True)
            bootstrap_distributions.append(bootstrap_sum)

        return np.array(bootstrap_distributions)

    total_deviances_per_sample = deviances.sum(axis = 0)

    # perform bootstrapping
    #bootstrap_distribution = bootstrap_deviances(deviances    = total_deviances_per_sample,
    #                                             n_bootstraps = n_bootstraps)

    def calculate_pvalues(total_deviance, bootstrap_dist):
        total_deviance = np.array(total_deviance).reshape(1,-1) # reshape to match dimenssions (1, n_samples)
        pvalues = np.mean(bootstrap_dist >= total_deviance,
                          axis = 0)
        return pvalues

    # compute pvalues for each sample
    pvalues = calculate_pvalues(total_deviance = total_deviances_per_sample,
                                bootstrap_dist = bootstrap_distribution)

    # create a summary DF to show each samples,total deviance and pvalue

    summary_df = pd.DataFrame({
        'Sample'        : data_df.columns,
        'Total_Deviance': total_deviances_per_sample,
        'Pvalue'        : pvalues})

    return summary_df

# For preprocessing count matrices to a genehist file
def CountsMatrixToGeneHist(df, binsize = .25, maxdepth = 18.5):

    # Remove rownames : TODO find numeric columns?
    df = df.iloc[:, 2:]
    df = df.apply(pd.to_numeric,
                  errors = 'coerce')
    # CPM transformation
    read_depth = df.sum(axis = 0) / 1e6
    df = df.div(read_depth,
                axis = 1)
    df = np.log2(df + 1)

    # Make bins 
    bins = np.arange(0,
                     maxdepth,
                     binsize)

    # Count how many genes are within the top and bottom of each bin
    bin_counts = []
    for bin in bins:
        top_bin = bin + binsize
        count_df = (df > bin) & (df <= top_bin)
        count_df_sum = count_df.sum(axis = 0)
        bin_counts.append(count_df_sum)

    # Transpose the matrix and add X-axis names
    bin_counts_df = pd.DataFrame(bin_counts)
    topbins = bins + binsize
    xaxis_names = [f"({bin},{topbin}]" for bin, topbin in zip(bins, topbins)]
    bin_counts_df.insert(loc    = 0,
                         column = 'Bins',
                         value  = xaxis_names)

    return bin_counts_df

# function designed to return pvalues for GCinformation if supplied
def GC_KSstats(_coverage_df):

    # initialize list to store values
    _kslst = []
    # calculate mean GBC for the whole library
    _mean_df = pd.DataFrame()
    _mean_df["gc_mean"] = _coverage_df.median(axis = 1)

    for column_name, _column_data in _coverage_df.iteritems():
        _ks_stat, _ks_pval = stats.ks_2samp(data1 = _column_data,
                                            data2 = _mean_df['gc_mean'])
        _kslst.append(_ks_stat)

    return _kslst

#  GeneBody Coverage Plot
def plotGC(SampleName,user_df, _coverage_df, _position,_figinfo,_fig=None):

    _axis = _fig.add_subplot(_figinfo["_subplot_rows"], 2, _position)
    ## Prep information for plotting
    BatchName = user_df.loc[user_df['Sample'] == SampleName,'Batch'].iloc[0]
    SamplesInBatch = user_df.loc[user_df['Batch'] == BatchName, 'Sample']

    # Calculate mean GeneBody Coverage for the entire library
    _mean_df = pd.DataFrame()
    _mean_df['gc_mean'] = _coverage_df[SamplesInBatch].median(axis=1)

    # acquire the pvalue information from the _figinfo object
    _ks_pval = user_df.loc[user_df['Sample'] == SampleName,'GBC_KSstats'].iloc[0]
    #Acquire the batch information from the user_df

    # Define the visual information for plotting each line
    sample_line_info = {'color'     : _figinfo["_curr_sample_color"],
                        'linewidth' : .5,
                        'linestyle' : '-'}
    ref_line_info = sample_line_info.copy()
    ref_line_info.update({'color' :  'lightgray',
                          'alpha'  : .4})
    mean_line_info = sample_line_info.copy()
    mean_line_info.update({'color'     : 'indigo',
                           'linestyle' : '--'})

    # Plot lines
    _x = np.arange(1, 101, 1) # The 1-100 GBC
    _axis.plot(_x,
               _coverage_df,
               **ref_line_info)
    _axis.plot(_x,
               _coverage_df[SampleName],
               **sample_line_info)
    _axis.plot(_x,
               _mean_df['gc_mean'],
               **mean_line_info)

    # include a 95% interval for each position
    _err = _coverage_df.std(axis=1)*2
    _axis.fill_between(_x,
                       _mean_df['gc_mean'] - _err,
                       _mean_df['gc_mean'] + _err,
                       facecolor = 'yellow',
                       alpha     =  0.3)

    # make the symbols for the legend
    _current_sample_line = matplotlib.lines.Line2D([0],
                                                   [0],
                                                   **sample_line_info)
    _background_lines = matplotlib.lines.Line2D([0],
                                                [0],
                                                **ref_line_info)
    _library_line = matplotlib.lines.Line2D([0],
                                            [0],
                                            **mean_line_info)

    _axis.legend([_current_sample_line,
                  _library_line,
                  _background_lines],
                 [f"Current Sample: K-S Stat = {_ks_pval:.2f}",
                  "Batch Mean",
                  "Batch Samples"],
                 loc      = 'lower center',
                 frameon  = False,
                 fontsize = _figinfo["_legend_size"],
                 ncol     = 1)
    _axis = mk_axes(_axis)
    _axis = needs_fail_or_warn(ax             = _axis,
                               current_sample = _ks_pval,
                               _figinfo       = _figinfo,
                               cutoff_key     = "GC_cutoff",
                               higher_lower   = "upper")

    _axis = set_ticks(_axis,
                      _figinfo["_tick_size"])
    _axis.set_xlim(0, 105)
    _axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_percent1))

    _axis.set_title("Gene Body Coverage",
                    fontsize = _figinfo["_title_size"])
    _axis.set_xlabel("Gene Body Position Percentile (5' " + u"\u2192" + " 3')",
                     fontsize = _figinfo["_label_size"])
    _axis.set_ylabel("Read Density",
                     fontsize = _figinfo["_label_size"])

    return _fig

def calcHistPval(_hist_df):

    # code for calculating Z value of number of expressed genes. In need of some improvement.
    _sum_df = _hist_df.sum().round()
    _zscore = stats.zscore(_sum_df)
    _pvals  = stats.norm.sf(abs(_zscore))

    return(_pvals)

# Not in use currently but in theory plotNegBin and plotGC could share some values although
#im not sure the added complexity is worth the squeeze.
class AbstractLinePlotter(ABC):

    def __init__(self,_ipTuple,_hist_df,_user_df,_position,_figinfo,_figure=None):
        self.IpTuple = _ipTuple
        self.Position= _position
        self.UserDf  = _user_df
        self.FigInfo = _figinfo
        self.Figure  = _figure
        self.HistDf  = _hist_df

    def AddLinePlot(self):
        _ax = self.Figure.add_subplot(self.FigInfo["_subplot_rows"],2,self.Position)

        _index_array = self.HistDf.iloc[:, 0]

        _low_vals =  []
        _high_vals = []

        for _i in _index_array:
            _low_vals.append(float(_i.strip('(').strip(']').split(',')[0]))
            _high_vals.append(float(_i.strip('(').strip(']').split(',')[1]))


# Plot 8 : Gene Expression Distribution Plot 
def plotNegBin(SampleName,UserDf, _hist_df, _position,_figinfo,_f=None):

    _ax = _f.add_subplot(_figinfo["_subplot_rows"],
                         2,
                         _position)
    _low_vals = []
    for _i in _hist_df.iloc[:, 0]:
        _low_vals.append(float(_i.strip('(').strip(']').split(',')[0]))

    ## Preparing the data_df and libMean_df for all bins
    _hist_df = _hist_df.drop(['Bins'],
                             axis = 1)
    _libMean_df = pd.DataFrame()

    # Get current batch
    current_batch = UserDf.loc[UserDf['Sample'] == SampleName,'Batch'].iloc[0]
    samples_in_current_batch = UserDf.loc[UserDf['Batch'] == current_batch, 'Sample'].tolist()
    batch_hist_df = _hist_df.loc[:, _hist_df.columns.isin(samples_in_current_batch)]


    # Get all samples that belong to current from UserDf
    _libMean_df['Mean'] = batch_hist_df.mean(numeric_only = True,
                                             axis         = 1)
    _current_samp_array = _hist_df[SampleName].values

    # code for calculating Z value of number of expressed genes. In need of some improvement.
    _sum_df = batch_hist_df.sum().round()
    _curr_sum = _current_samp_array.sum().round()
    _curr_ndx = np.where(_sum_df == _curr_sum)[0][0]
    hist_mean = _sum_df.mean().round()
    _zscore = stats.zscore(_sum_df)
    _pvals  = stats.norm.sf(abs(_zscore))
    _curr_pval = _pvals[_curr_ndx]

    # Define the visual information for plotting each line
    sample_line_info = {'color'     : _figinfo["_curr_sample_color"],
                        'linewidth' : .5,
                        'linestyle' : '-'}
    ref_line_info = sample_line_info.copy()
    ref_line_info.update({'color' :  'lightgray',
                          'alpha'  : .4})
    mean_line_info = sample_line_info.copy()
    mean_line_info.update({'color'     : 'indigo',
                           'linestyle' : '--'})

    _ax.plot(_low_vals,
             _hist_df,
             **ref_line_info)
    _ax.plot(_low_vals,
             _hist_df[SampleName],
             zorder = 24,
             **sample_line_info)
    _ax.plot(_low_vals,
             _libMean_df['Mean'],
             zorder    = 23,
             **mean_line_info)


    _current_samp_line = matplotlib.lines.Line2D([0],
                                                 [0],
                                                 **sample_line_info)
    _lib_line = matplotlib.lines.Line2D([0],
                                        [0],
                                        **mean_line_info)

    _ax.legend([_current_samp_line,
                _lib_line],
               [f'Current Sample: {_curr_sum.item()} Detected Genes (p={_curr_pval.item():.3f})',
                f'Batch Mean: {hist_mean} Detected Genes'],
               loc      = 'upper right',
               frameon  = False,
               fontsize = _figinfo["_legend_size"],
               ncol     = 1)

    _ax = set_ticks(_ax,_figinfo["_tick_size"])

    _ax.set_xlim(0, 10)

    _ax.set_title("Gene Expression Distribution",
                  fontsize=_figinfo["_title_size"])
    _ax.set_xlabel("Expression Level (log2(CPM)+1)",
                   fontsize=_figinfo["_label_size"])
    _ax.set_ylabel("# Genes",
                   fontsize=_figinfo["_label_size"] )

    _ax = mk_axes(_ax)
    _ax = needs_fail_or_warn(ax             = _ax,
                             current_sample = _curr_pval,
                             _figinfo       = _figinfo,
                             cutoff_key     = "_alpha",
                             higher_lower   = "lower")

    return _f

def get_axis_range(axis_lim):
     """
     Get the range size  of an axis matplotlib lim object
     Input:
         axis_lim: a matplotlib lim object
     Output: The size of the range
     """
     rng = axis_lim[1] - axis_lim[0]
     return rng

def repl_missing_values_indict(_indict,_repldict):
    for key, value in _indict.items():
    # If the value is NaN, replace it with the corresponding automatically generated keys
        if value != value:
            _indict.update({key: _repldict[key]})

def read_file(filepath, **kwargs):
    """
    Reads a file and returns a DataFrame.
    Supports both CSV and Excel file formats.

    Parameters:
    - filepath (str): Path to the file (CSV or Excel).
    - kwargs: Additional arguments passed to the respective pandas read functions.

    Returns:
    - pd.DataFrame: DataFrame containing the file's data.
    """
    if filepath.endswith(".csv"):
        return pd.read_csv(filepath, **kwargs)
    elif filepath.endswith((".xls", ".xlsx")):
        return pd.read_excel(filepath, **kwargs)
    else:
        raise ValueError("Unsupported file format. Please provide a CSV or Excel file.")

if __name__ == "__main__":
    main()
