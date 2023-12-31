import matplotlib
from abc import ABC, abstractmethod
import matplotlib.gridspec as gridspec
import os
import argparse
matplotlib.use('PDF')
import seaborn
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
from statsmodels.stats.weightstats import ztest
import matplotlib.patches as mpatches
from Sam_PyUtils import *
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
import bootstrapped.compare_functions as bs_compare
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
    # store the library of methods that different subclasses use depending on their situation.   
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
        print(self.mapping_names)
        print(self.input_df.columns)
        self.input_df.rename(columns=self.mapping_names, inplace=True)
        print(self.input_df.columns)
    def _transform_values(self):
        """ Calculate columns that are transformations of the supplied inputs """
        self.input_df["% Uniquely Aligned Reads"] = (
        (self.input_df["# Uniquely Aligned Reads"]) / (self.input_df["Sequencing Depth"] * self.input_df["% of Reads After Trimming"])) * 10000

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
        "Sequencing Depth",
        "% of Reads After Trimming",
        "# Uniquely Aligned Reads",
        "% Uniquely Aligned Reads",
        "% Exonic Reads / Aligned Reads",
        "# rRNA Reads",
        "% Overrepresented Sequences (Pre-trim)",
        "% Adapter Content (Pre-trim)",
        "% Overrepresented Sequences (Post-trim)",
        "% Adapter Content (Post-trim)",
        "Project",
        "Batch"]
    input_program_names = [
        "Sample",
        "Input_Size",
        "Percent_PostTrim",
        "Num_Uniquely_Aligned",
        "Percent_Uniquely_Aligned",
        "Percent_Exonic",
        "Num_Uniquely_Aligned_rRNA",
        "Percent_Overrepresented_Seq_Untrimmed",
        "Percent_Adapter_Content_Untrimmed",
        "Percent_Overrepresented_Seq_Trimmed",
        "Percent_Adapter_Content_Trimmed",
        "Project",
        "Batch"]

    def adapt_input(self):
        self.make_mapping_names()
        self._validate_column_names()
        self._fill_missing_values()
        self._transform_values()
        self._change_column_names()
        self._reorder_columns()       

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
        "% Overrep Sequences (Pre-Trim)",
        "% Adapter Content (Pre-Trim)",
        "% Overrep Sequences (Post-Trim)",
        "% Adapter Content (Post-Trim)",
        "# Detected Genes",
        "Gene Body Coverage (Pval)"]

    input_program_names = [
        "cutoff",
        "_ipReads_cutoff",
        "_trimmedReads_cutoff",
        "_uniqAligned_cutoff",
        "_exonMapping_cutoff",
        "_riboScatter_cutoff",
        "_violin_cutoff_overrep_untrimmed",
        "_violin_cutoff_adapter_untrimmed",
        "_violin_cutoff_overrep_trimmed",
        "_violin_cutoff_adapter_trimmed",
        "Dist_of_gene_expression",
        "GeneBody_Coverage"]

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

        ax.plot(cutoff,ax.get_ylim()[1] - (18 * plot_scalar), marker='v', ms=0.8, c=clr)
        ax.text(cutoff,ax.get_ylim()[1] - (13 * plot_scalar), txt, fontsize=4, color=clr,
                 horizontalalignment='center')                                                                  
        return ax 
                                                                                            
    plot_scalar = get_axis_range(ax.get_ylim()) / 100

    ax = add_vert_marker(ax,_figinfo["_fail_cutoffs"][cutoff_key],plot_scalar,_figinfo["_fail_color"],"Fail")
    ax = add_vert_marker(ax,_figinfo["_warn_cutoffs"][cutoff_key],plot_scalar,_figinfo["_warn_color"],"Warn")
    
    return ax
   
# This function calculates whether the current sample label orientation needs 
# to be adjusted and returns the required vars
def adjust_flag(_ax,_current_sample,_lib_mean):
    if _current_sample >=  _lib_mean:
        # plot the current sample line
        _ax.text(_current_sample + (get_axis_range(_ax.get_xlim())/ 100), (_ax.get_ylim()[1] / 2), '{:2.2f}M'.format(_current_sample),
                 rotation= 270, fontsize=3, zorder=2)
        # plot the library mean line
        _ax.text(_lib_mean - (get_axis_range(_ax.get_xlim())/30),((_ax.get_ylim()[1] / 2) + 1), '{:2.2f}M'.format(_lib_mean), rotation=90,
                 fontsize=3, zorder=2)
    else:
        # plot the current sample line
        _ax.text(_current_sample - (get_axis_range(_ax.get_xlim())/50),(_ax.get_ylim()[1] / 2), '{:2.2f}M'.format(_current_sample),
                 rotation= 90, fontsize=3, zorder=2)
        # plot the library mean line
        _ax.text(_lib_mean + (get_axis_range(_ax.get_xlim())/80), ((_ax.get_ylim()[1] / 2) + 1), '{:2.2f}M'.format(_lib_mean), rotation= 270,
                 fontsize=3, zorder=2)
    return _ax

# set the legend for figures with warn / fail IE 1_6)
def legend_setup_1_6(_ax,_line1,_line2,_figinfo,cutoff_key,_loc):
    _fail_label = mpatches.Patch(color=_figinfo["_fail_color"], label='Fail Cutoff')   
    _warn_label=  mpatches.Patch(color=_figinfo["_warn_color"], label='Warn Cutoff')    
    _ax.legend([_line1,
                    _line2,
                    _fail_label,
                    _warn_label],
                ["Current Sample", 
                    "Batch Mean",
                    "Fail (" + str(round(_figinfo["_fail_cutoffs"][cutoff_key],2)) + ")",
                    "Warn (" + str(round(_figinfo["_warn_cutoffs"][cutoff_key],2)) + ")"],
                loc= _loc,
                frameon=False,
                fontsize=_figinfo["_legend_size"])
    return _ax 

# the goal of this function is to set which axes of the regular axis and the Kernel density axis 
def mk_axes(_plt_ax,_kd_ax = None):
    
    for label in (_plt_ax.get_xticklabels() + _plt_ax.get_yticklabels()):
        label.set_fontsize(4)
    # set the plot axis
    _plt_ax.spines['top'].set_visible(False)
    _plt_ax.spines['right'].set_visible(False)
    _plt_ax.spines['left'].set_visible(True)
    _plt_ax.spines['bottom'].set_visible(True)
    _plt_ax.spines['left'].set_color('black')
    _plt_ax.spines['bottom'].set_color('black')

    _plt_ax.spines['left'].set_linewidth(0.55)
    _plt_ax.spines['bottom'].set_linewidth(0.55)
    _plt_ax.set_facecolor('white')
    
    if _kd_ax != None: 
        # set the kernel density axis to be invisible 
        _kd_ax.yaxis.set_ticks([])
        _kd_ax.yaxis.label.set_visible(False)

        _kd_ax.spines['top'].set_visible(False)
        _kd_ax.spines['right'].set_visible(False)
        _kd_ax.spines['bottom'].set_visible(False)
        _kd_ax.spines['left'].set_visible(False)

        return _plt_ax,_kd_ax
    else:
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
                font_size = 4
                background_color = 'tomato'
                text_color = 'yellow'
                font_stretch = 'extra-expanded'
            elif flag_type == "warn":
                text = "WARNING"
                font_size = 3.7
                background_color = 'yellow'
                text_color = 'red'
                font_stretch = 'expanded'
            else:
                raise ValueError("Invalid flag type")

            anch_text = matplotlib.offsetbox.AnchoredText(text, pad=0.0001, borderpad=1,
                                                          loc=3,
                                                          prop=dict(size=font_size, snap=True, backgroundcolor=background_color,
                                                                    color=text_color, alpha=0.9, zorder=5, fontweight='roman',
                                                                    fontfamily='serif', fontstretch=font_stretch),
                                                          frameon=True,
                                                          bbox_to_anchor=(0., 1.03),
                                                          bbox_transform=self.ax.transAxes)

            self.ax.add_artist(anch_text)
            return None
 
    def insert_flag(ax, cutoff, flag_func):
        if higher_lower == "lower" and current_sample <= cutoff  or higher_lower == "upper" and current_sample >= cutoff:
            flag_func(ax)

    flag_inserter = FlagInserter()
    insert_flag(ax, flag_inserter.cutoff_warn, lambda ax: flag_inserter.make_flag("warn"))
    insert_flag(ax, flag_inserter.cutoff_fail, lambda ax: flag_inserter.make_flag("fail"))
    
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

# The objective of this function is to generate the dynamic fail cutoffs for sample display, based on background data

class CutoffCalculator:
    def __init__(self, bgd_df, alph):
        self.bgd_df = bgd_df
        self.onesided_alph = 2*alph
        self.cutoffs_dict = {}
        self.cutoffs_dict["_alpha"] = alph

    def __call__(self):
        self.calculate_cutoff("Input_Size", "lower", "_ipReads_cutoff")
        self.calculate_cutoff("Percent_PostTrim", "lower", "_trimmedReads_cutoff")
        self.calculate_cutoff("Percent_Uniquely_Aligned", "lower", "_uniqAligned_cutoff")
        self.calculate_cutoff("Percent_Exonic", "lower", "_exonMapping_cutoff")
        self.calculate_cutoff_ratio("Num_Uniquely_Aligned_rRNA", "Num_Uniquely_Aligned", "upper", "_riboScatter_cutoff")
        self.calculate_cutoff("Percent_Overrepresented_Seq_Untrimmed", "upper", "_violin_cutoff_overrep_untrimmed")
        self.calculate_cutoff("Percent_Adapter_Content_Untrimmed", "upper", "_violin_cutoff_adapter_untrimmed")
        self.calculate_cutoff("Percent_Overrepresented_Seq_Trimmed", "upper", "_violin_cutoff_overrep_trimmed")
        self.calculate_cutoff("Percent_Adapter_Content_Trimmed", "upper", "_violin_cutoff_adapter_trimmed")
        
        return self.cutoffs_dict

    def calculate_cutoff_ratio(self, column1, column2, upper_lower, cutoff_name):
        vec = np.array(self.bgd_df.loc[:, column1] / self.bgd_df.loc[:, column2])
        bootstrap_mean = bs.bootstrap(vec, stat_func=bs_stats.mean).value
        bootstrap_std  = bs.bootstrap(vec, stat_func=bs_stats.std).value
        conf_size = norm.ppf(self.onesided_alph) * bootstrap_std
        if upper_lower == "upper":
            self.cutoffs_dict[cutoff_name] = bootstrap_mean - conf_size
        if upper_lower == "lower":
            self.cutoffs_dict[cutoff_name] = bootstrap_mean + conf_size

    def calculate_cutoff(self,column,upper_lower,cutoff_name):
        vec = np.array(self.bgd_df.loc[:, column])
        bootstrap_mean = bs.bootstrap(vec, stat_func=bs_stats.mean).value
        bootstrap_std  = bs.bootstrap(vec, stat_func=bs_stats.std).value
        conf_size = norm.ppf(self.onesided_alph) * bootstrap_std
        if upper_lower == "upper":
            self.cutoffs_dict[cutoff_name] = bootstrap_mean - conf_size
        if upper_lower == "lower":
            self.cutoffs_dict[cutoff_name] = bootstrap_mean + conf_size


def gen_cutoffs(bgd_df, alph):
    calculator = CutoffCalculator(bgd_df, alph)
    return calculator()


def fmt_scatter_million(_x, _pos):
    return '%1.1f' % (_x * 1e-6)

def fmt_million(_x, _pos):
    return '{0:.0f}'.format(_x)

def fmt(_x, _pos):
    return '{0:.1f}'.format(_x)

def values_to_percentiles(values):
    """
    Convert a vector of values to their associated percentile ranks on a standard normal distribution.
    
    Parameters
    ----------
    values : list or numpy array of float
        A vector of values you want to convert to their associated percentile ranks.
    """
    
    # Convert the input values to a numpy array if not already
    values = np.asarray(values)

    # Calculate the mean and standard deviation of the input values
    mean = np.mean(values)
    std_dev = np.std(values)

    # Standardize the values by subtracting the mean and dividing by the standard deviation
    standardized_values = (values - mean) / std_dev

    # Calculate the percentile rank for each standardized value using the cumulative distribution function (CDF)
    percentiles = norm.cdf(standardized_values)

    return percentiles

def set_ticks(_ax,_tick_size):
    _ax.tick_params(axis='x', which='both', length=1, width=0.5, labelbottom=True, bottom=True, labelsize= _tick_size,
                      direction='out', pad=2)
    _ax.tick_params(axis='y', which='both', length=1, width=0.5, labelsize= _tick_size, labelleft=True, left=True,
                      direction='out', pad=2)   

    return _ax

def label_anno(ax, line, label, color='0.5', fs=3, halign='left', valign='center_baseline'):

    xdata, ydata = line.get_data()
    x1 = xdata[0]
    x2 = xdata[-1]
    y1 = ydata[0]
    y2 = ydata[-1]

    if halign.startswith('l'):
        xx = x1
        halign = 'left'
    elif halign.startswith('r'):
        xx = x2
        halign = 'right'
    elif halign.startswith('c'):

        if ax.get_xscale() == 'log':
            xx = 10 ** (0.5 * (np.log10(x1) + np.log10(x2)))
        else:
            xx = 0.5 * (x1 + x2)
        halign = 'center'
    else:
        raise ValueError("Unrecognized `halign` = '{}'.".format(halign))

    if ax.get_xscale() == 'log' and ax.get_yscale() == 'log':
        yy = 10 ** (np.interp(np.log10(xx), np.log10(xdata), np.log10(ydata)))
    elif ax.get_xscale() == 'log' and ax.get_yscale() != 'log':
        yy = np.interp(np.log10(xx), np.log10(xdata), ydata)
    elif valign.startswith('t'):
        valign = 'top'
        yy = np.interp(xx, xdata, ydata)
    else:
        yy = np.interp(xx, xdata, ydata)


    ylim = ax.get_ylim()
    xytext = (0, 0)
    text = ax.annotate(label, xy=(xx, yy), xytext=xytext, textcoords='offset points', size=fs, color=color, zorder=1, horizontalalignment=halign, verticalalignment=valign)

    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))

    rise = (sp2[1] - sp1[1])
    run = (sp2[0] - sp1[0])

    slope_degrees = np.degrees(np.arctan2(rise, run))
    text.set_rotation_mode('anchor')
    text.set_rotation(slope_degrees)
    ax.set_ylim(ylim)

    return text

# takes two vectors, finds the maximum and minimum, and creates the number of bins specified
def make_bins(vec1,vec2,num_bins):

    plt_min,plt_max = pd.concat([vec1,vec2]).agg(['min','max'])
    bins = np.linspace(plt_min,plt_max, num_bins + 1)
    
    return bins

def mkTitlePage(_figinfo):

    """ 
    Make QC heatmap data
    The Goal of this function is to calculate a matrix that can be used with matplotlibs heatmap
    function which contains whether a sample passed, failed, or was warned for each test
    """
    
    
    def mk_cutoff_descript(descriptor,cutoff_set):
        
        cutoff_descriptor = (descriptor +': Default alpha =' + cutoff_set["_alpha"]  +
                                '| Sequencing Depth = ' + '{:.3f}'.format(cutoff_set["_ipReads_cutoff"],3)) +
                                '| Trimming = ' + '{:.3f}'.format(cutoff_set["_trimmedReads_cutoff"],3)) +  
                                '| Alignment = ' + '{:.3f}'.format(cutoff_set["_uniqAligned_cutoff"],3)) +  
                                '| Gene Exon Mapping = ' + '{:.3f}'.format(cutoff_set["_exonMapping_cutoff"],3)) + 
                                '\n| Ribosomal RNA = ' + '{:.3f}'.format(cutoff_set["_riboScatter_cutoff"],3)) + 
                                '| Adapter Contamination = ' + '{:.3f}'.format(cutoff_set["_violin_cutoff_adapter_trimmed"],3)) + 
                                '| Overrep. Seq  Contamination = ' + '{:.3f}'.format(cutoff_set["_violin_cutoff_overrep_trimmed"],3)) + 
                                '| Gene Body Coverage = ' + cutoff_set["_alpha"] + 
                                '| Detected Genes = ' + cutoff_set["_numGene_cutoff"])
        return(cutoff_descriptor)  

    fig = plt.figure()

    # add text
    fig.text(.5,.965,"QC Plotter Input Summary",ha='center',va='top',fontsize = 14)
    fig.text(.005,.985,datetime.now().strftime("%d/%m/%Y %H:%M:%S"),fontsize = 4)
    fig.text(.005,.97,"Version 2.0",fontsize = 4)
    fig.text(.5,.5,"Input table =     " + _figinfo["_ip_filename"] + "\nOutput location =    " + \
    _figinfo["_op_filename"] + "\nBackground table =    " + _figinfo["_bgd_filename"] + 
    "\nGene Body Coverage file =    " + str(_figinfo["_gc_file"]) + \
    "\nGene read depth distribution histogram file =" + \
    str(_figinfo["_hist_file"]) + "\nCutoff file =    " + str(_figinfo["_cutoff_filename"]),
    fontsize = 6,ha="center",va="center")

    # show cutoffs

    warn_descript = mk_cutoff_descript("Warn cutoffs",_figinfo["_warn_cutoffs"]) 
    fail_descript = mk_cutoff_descript("Fail cutoffs",_figinfo["_fail_cutoffs"]) 


    fig.text(s= warn_descript , 
                    x = .5, y = .1, fontsize = 5,
                     ha='center', va='top',fontweight='book', style = 'italic')


    fig.text(s= fail_descript , 
                    x = .5, y = .05, fontsize = 5,
                     ha='center', va='top',fontweight='book', style = 'italic')

    return fig


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
        ("Input_Size", "_ipReads_cutoff", lower_status_strategy),
        ("Percent_PostTrim", "_trimmedReads_cutoff", lower_status_strategy),
        ("Percent_Uniquely_Aligned","_uniqAligned_cutoff",lower_status_strategy),
        ("Percent_Exonic","_riboScatter_cutoff",lower_status_strategy),
        ("Num_Uniquely_Aligned_rRNA","_riboScatter_cutoff",upper_status_strategy), # HOLD UNTIL I FIX THE ISSUES WITH RIBO Pretty sure this is fixed actually
        ("Percent_Overrepresented_Seq_Trimmed","_violin_cutoff_overrep_trimmed",upper_status_strategy),
        ("Percent_Adapter_Content_Trimmed","_violin_cutoff_adapter_trimmed",upper_status_strategy)]
    if _figinfo["_hist_exists"]:
        strategies.append(("_hist_pvals","_alpha",lower_status_strategy))
    if _figinfo["_gbc_exists"]:
        strategies.append(("_gbc_pvals","_alpha",lower_status_strategy))    

    _htmat = np.zeros((len(_userDf), 9))

    for _tuple in _userDf.itertuples():
        for i, (column, key, strategy) in enumerate(strategies):
        
            if column == "Num_Uniquely_Aligned_rRNA":
                test_value = _userDf.iloc[_tuple.Index][column] / _userDf.iloc[_tuple.Index]["Num_Uniquely_Aligned"]
                print("rRNA test value is",test_value)
            else:
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
    colors = ["grey","goldenrod","red"]
    cm = matplotlib.colors.ListedColormap(colors)
    sample_names = heatmap_data.Sample
    heatmap_data = heatmap_data.drop("Sample",axis = 1)
    fig2,ax = plt.subplots()
    fig2.text(s= "Summary of QC Metrics",x = .5,y = .9,fontsize = 10,ha = 'center')  
    seaborn.heatmap(heatmap_data.values,ax=ax,
                    xticklabels=["Sequencing Depth","Trimming","Alignment","Exon Mapping","Ribosomal RNA",
                                  "Sequence Contamination (Overrep)","Sequence Contamination (Adapter)",
                                 "# Detected Genes","Gene Body Coverage"],
                    yticklabels=sample_names,
                    cmap = cm)
    # change y-axis tick label font size
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

    plt.subplots_adjust(left=0.3, bottom=0.3, right=0.7, top=0.8)
    
    return fig2



def plotHist_ipSize(_ip_tuple, _userDf, _background_df, _position,_figinfo,_figure=None):

    axis = _figure.add_subplot(_figinfo["_subplot_rows"], 2, _position)
 
    _xmin,_xmax = _background_df.loc[:,"Input_Size"].agg(["min","max"])
    _bins = make_bins(_background_df.loc[:,"Input_Size"],_user_df.loc[:,"Input_Size"],_figinfo["_bin_num"])
    
    _lib_mean = _user_df.loc[:, 'Input_Size'].mean()
    _current_sample = _user_df.loc[:,'Input_Size'][_ip_tuple.Index]

    _background_df["Input_Size"].plot(kind='hist', bins=_bins, ax=axis, color='lightgray')

    axis1 = axis.twinx()
    sns.distplot(_background_df["Input_Size"], hist=False, bins=_bins, ax=axis1, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})
    axis.set_xlim(_xmin,_xmax)
    axis = set_ticks(axis,_figinfo["_tick_size"])

    axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_million))

    axis.set_title("Sequencing Depth",fontsize = _figinfo["_title_size"])
    axis.set_xlabel('Total Reads (Millions)', labelpad=1, fontsize= _figinfo["_label_size"])

    axis.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    axis.set_ylabel('Frequency', labelpad=2, fontsize= _figinfo["_label_size"])

    axis = adjust_flag(axis,_current_sample,_lib_mean)

    ### Adding cutoff markers
    axis = add_warn_fail_markers(_figinfo,axis,"_ipReads_cutoff")

    # Current Sample Line and Label
    _line1 = axis.axvline(x=_current_sample, alpha=0.8, color=_figinfo["_curr_sample_color"], linestyle='-', linewidth=0.5,
                         label='{:2.2f}M'.format(_current_sample))

    # Current Library Mean Line and Label
    _line2 = axis.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5,
                         label='{:2.2f}M'.format(_lib_mean))

    _kde_line = matplotlib.lines.Line2D([0], [0], color="gray", linewidth=0.5, linestyle='-')

    # set up axes
    axis = legend_setup_1_6(axis,_line1,_line2,_figinfo,"_ipReads_cutoff","upper left")

    #set axes to be visible or not
    axis,axis1 =  mk_axes(axis,axis1)
    axis = needs_fail_or_warn(axis,_current_sample,_figinfo,"_ipReads_cutoff","lower")

    return _figure


# Plot 2 : Trimming Percentage
def plotHist_trimming(_ip_tuple, _user_df, _background_df, _colname, _position,_figinfo,_figure=None):
  
    axis = _figure.add_subplot(_figinfo["_subplot_rows"],2, _position)    

    _xmin,_xmax = _background_df.loc[:,_colname].agg(["min","max"])

    _bins = make_bins(_background_df.loc[:,_colname],_user_df.loc[:,_colname],_figinfo["_bin_num"])

    _current_sample = _user_df.loc[:,_colname][_ip_tuple.Index]
    _lib_mean = _user_df[_colname].mean()

    axis.hist(x=_background_df[_colname], bins=_bins, histtype='bar', color='lightgray')

    axis1 = axis.twinx()
    sns.distplot(_background_df[_colname], hist=False, bins=_bins, ax=axis1, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})
    
    axis.set_xlim(_xmin, _xmax)
    axis = set_ticks(axis,_figinfo["_tick_size"])
    axis.set_title("Trimming", fontsize= _figinfo["_title_size"])

    axis.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    axis.set_ylabel('Frequency', labelpad=2, fontsize= _figinfo["_label_size"])
    axis.set_xlabel('% Post-Trim / Total Reads', labelpad=1, fontsize= _figinfo["_label_size"])

    ### Adding cutoff markers
    axis = add_warn_fail_markers(_figinfo,axis,"_trimmedReads_cutoff")
    
    axis = adjust_flag(axis,_current_sample,_lib_mean)

    # Current Sample Line and Label
    _line1 = axis.axvline(x=_current_sample, alpha=0.8, color=_figinfo["_curr_sample_color"], linestyle='-', linewidth=0.5,
                           label='{:2.2f}%'.format(_current_sample))

    # Current Library Mean Line and Label
    _line2 = axis.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5,
                           label='{:2.2f}%'.format(_lib_mean))

    ## Superimpose the Kernel Density Estimate line over the distribution
    _kde_line = matplotlib.lines.Line2D([0], [0], color="dimgray", linewidth=0.5, linestyle='-')

    axis = legend_setup_1_6(axis,_line1,_line2,_figinfo,"_trimmedReads_cutoff","upper left")

    axis,axis1 =  mk_axes(axis,axis1)
    axis = needs_fail_or_warn(axis,_current_sample,_figinfo,"_trimmedReads_cutoff","lower")

    return _figure


#### Plot 3: Alignment Percentage ####
def plotHist_alignment(_ip_tuple, _user_df, _background_df, _colname, _plot_label,_figinfo,_figure=None):
   
    _bins = np.arange(0, 100 + 1,101 / _figinfo["_bin_num"] )

    _current_sample = _user_df.loc[:,_colname][_ip_tuple.Index]

    _lib_mean = _user_df[_colname].mean()
    _axis_plt3 = _figure.add_subplot(_figinfo["_subplot_rows"], 2, _position)

    _axis_plt3.hist(x=_background_df[_colname], bins=_bins, histtype='bar', color='lightgray')

    _axis1_plt3 = _axis_plt3.twinx()
    sns.distplot(_background_df[_colname], hist=False, bins=_bins, ax=_axis1_plt3, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})


    _axis_plt3.set_xlim(0, 105)

    _axis_plt3 = set_ticks(_axis_plt3,_figinfo["_tick_size"])

    _axis_plt3.set_title("Alignment", fontsize=_figinfo["_title_size"])

    # Set labels
    _axis_plt3.set_xlabel('% Uniquely Aligned / Post-Trim Reads', labelpad=1, fontsize= _figinfo["_label_size"])
    _axis_plt3.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    _axis_plt3.set_ylabel('Frequency', labelpad=2, fontsize= _figinfo["_label_size"])

    ### Adding cutoff markers
    _axis_plt3 = add_warn_fail_markers(_figinfo,_axis_plt3,"_uniqAligned_cutoff")
    _axis_plt3 = adjust_flag(_axis_plt3,_current_sample,_lib_mean)

    # Current Sample Line and Label
    _line1 = _axis_plt3.axvline(x=_current_sample, alpha=0.8, color=_figinfo["_curr_sample_color"], linestyle='-', linewidth=0.5, label='{:2.2f}%'.format(_current_sample))

    # Current Library Mean Line and Label
    _line2 = _axis_plt3.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5, label='{:2.2f}%'.format(_lib_mean))

    ## Superimpose the Kernel Density Estimate line over the distribution
    _kde_line = matplotlib.lines.Line2D([0], [0], color="dimgray", linewidth=0.5, linestyle='-')

    _axis_plt3 = legend_setup_1_6(_axis_plt3,_line1,_line2,_figinfo,"_uniqAligned_cutoff","upper left")

    _axis_plt3,_axis1_plt3 =  mk_axes(_axis_plt3,_axis1_plt3)
    _axis_plt3 = needs_fail_or_warn(_axis_plt3,_current_sample,_figinfo,"_uniqAligned_cutoff","lower")
    
    return _figure



#### Plot 4: Gene Exon Mapping ####
def plotHist_exonMapping(_ip_tuple, _user_df, _background_df, _colname,_position,_figinfo,_figure=None):
    
    _bins = np.arange(0, 100 + 1, 101/_figinfo["_bin_num"])

    _user_minusBatchMean_df = _user_df.drop(_user_df.tail(1).index)

    _current_sample = _ip_tuple.Percent_Exonic
    _lib_mean = _user_minusBatchMean_df[_colname].mean()
    
    _axis_plt4 = _figure.add_subplot(_figinfo["_subplot_rows"], 2, _position)

    _axis_plt4.hist(x=_background_df[_colname], bins=_bins, histtype='bar', color='lightgray')

    _axis1_plt4 = _axis_plt4.twinx()
    sns.distplot(_background_df[_colname], hist=False, bins=_bins, ax=_axis1_plt4, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})


    _axis_plt4.set_xlim(0, 105)

    _axis_plt4 = set_ticks(_axis_plt4,_figinfo["_tick_size"])
    _axis_plt4.set_title("Exon Mapping", fontsize=_figinfo["_title_size"])

    _axis_plt4.set_xlabel('% Mapped / Aligned Reads', labelpad=1, fontsize=_figinfo["_label_size"])
    _axis_plt4.set_ylabel('Frequency', labelpad=2, fontsize= _figinfo["_label_size"])

    ### Adding cutoff markers
    _axis_plt4 = add_warn_fail_markers(_figinfo,_axis_plt4,"_exonMapping_cutoff")

    _axis_plt4.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    _axis_plt4 = adjust_flag(_axis_plt4,_current_sample,_lib_mean)

    # Current Sample Line and Label
    _line1 = _axis_plt4.axvline(x=_current_sample, alpha=0.8, color=_figinfo["_curr_sample_color"], linestyle='-', linewidth=0.5, label='{:2.2f}%'.format(_current_sample))

    # Current Library Mean Line and Label
    _line2 = _axis_plt4.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5, label='{:2.2f}%'.format(_lib_mean))

    ## Superimpose the Kernel Density Estimate line over the distribution
    _kde_line = matplotlib.lines.Line2D([0], [0], color="dimgray", linewidth=0.5, linestyle='-')


    _axis_plt4 = legend_setup_1_6(_axis_plt4,_line1,_line2,_figinfo,"_exonMapping_cutoff","upper left")
    _axis_plt4,_axis1_plt4 =  mk_axes(_axis_plt4,_axis1_plt4)
    _axis_plt4 = needs_fail_or_warn(_axis_plt4,_current_sample,_figinfo,"_exonMapping_cutoff","lower")
    
    return _figure

#### Plot 5: rRNA Scatter ####
def plotScatter_rRNA(_in_tup, _userDf, _background_df, _pos,_figinfo,_f=None):
    _plotter_df = pd.concat([_background_df, _userDf],sort = True)

    # Assign color for current project's library (all samples in the current project)
    _plotter_df["scatter_color"] = np.where(_plotter_df["Sample"].isin(_userDf["Sample"]), "indigo", "lightgray")

    # Assign separate color for current sample on each page
    _plotter_df.loc[_plotter_df["Sample"] == _in_tup[1], "scatter_color"] = _figinfo["_curr_sample_color"]

    _ax = plt.subplot(_figinfo["_subplot_rows"], 2, _pos)


    ## Regression line (gradient slope)
    X = _plotter_df.loc[:, "Num_Uniquely_Aligned"].values.reshape(-1, 1)
    Y = _plotter_df.loc[:, "Num_Uniquely_Aligned_rRNA"].values.reshape(-1, 1)
    linear_regressor = LinearRegression()
    linear_regressor.fit(X, Y)
    Y_pred = linear_regressor.predict(X)

    _ax.plot(X, Y_pred, c='black', linewidth=0.7, linestyle='-', alpha=1)
    _ax.scatter(x=_plotter_df['Num_Uniquely_Aligned'], y=_plotter_df['Num_Uniquely_Aligned_rRNA'], s=0.8, c=_plotter_df["scatter_color"])

    #separate scatter call for the sample so it can have a nique size and shape
    _intupdf = _plotter_df.loc[_plotter_df["Sample"] == _in_tup[1]]
    _ax.scatter(x=_intupdf['Num_Uniquely_Aligned'], y=_intupdf['Num_Uniquely_Aligned_rRNA'],marker = "*", s=20, c=_intupdf["scatter_color"])
    _ax.set_title("Ribosomal RNA", fontsize= _figinfo["_title_size"])
    _ax = set_ticks(_ax,_figinfo["_tick_size"])

    _ax.set_xlabel("Total Uniquely Aligned Reads (Millions)", fontsize= _figinfo["_label_size"], labelpad=2)
    _ax.set_ylabel("Aligned rRNA Reads (Millions)", fontsize= _figinfo["_label_size"], labelpad=2)

    x_bottom, x_top = plt.xlim()
    y_bottom, y_top = plt.ylim()

    # Plotting the ratio line 
    _slope_current = float(_in_tup[7] / _in_tup[4])

    xmin, xmax = _ax.get_xlim()
    line_x0 = 0
    line_x1 = xmax

    line_y0 = 0
    line_y1_warn = _figinfo["_warn_cutoffs"]["_riboScatter_cutoff"] * (line_x1 - line_x0) + line_y0
    line_y1_fail = _figinfo["_fail_cutoffs"]["_riboScatter_cutoff"] * (line_x1 - line_x0) + line_y0

    _ax.plot([line_x0, line_x1], [line_y0, line_y1_warn], c=_figinfo["_warn_color"], linewidth=1, linestyle='--', alpha=0.3, label="Warn")
    _ax.plot([line_x0, line_x1], [line_y0, line_y1_fail], c=_figinfo["_fail_color"], linewidth=1, linestyle='--', alpha=0.3, label="Fail")

    # Set axes margins for padding on both axes
    _ax.margins(0.01)

    _ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2500000))
    _ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_scatter_million))

    #_ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2500000))
    _ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_scatter_million))

    _ax.set_aspect('auto', adjustable='box', anchor='SW')

    _historic_data = matplotlib.lines.Line2D([0], [0], color='w', markerfacecolor='darkgray', marker='o', linewidth=1, markersize=3.5)
    _curr_lib = matplotlib.lines.Line2D([0], [0], color='w', markerfacecolor='indigo', marker='o', linewidth=1, markersize=3.5)
    _regression_gradient = matplotlib.lines.Line2D([0], [0], color='black', linewidth=0.6)
    _curr_samp = matplotlib.lines.Line2D([0], [0], color='w', markerfacecolor=_figinfo["_curr_sample_color"], marker='*', linewidth=1, markersize=6)
    _mean_label = mpatches.Patch(color='black', label='Mean Slope')   
    _fail_label = mpatches.Patch(color=_figinfo["_fail_color"], label='Fail Cutoff')   
    _warn_label=  mpatches.Patch(color=_figinfo["_warn_color"], label='Warn Cutoff')    
    _ax.legend([_curr_samp,
                    _curr_lib,
                    _fail_label,
                    _warn_label,
                    _mean_label],
                ["Current Sample", 
                    "Batch Samples",
                    "Fail (" + "{:.0%}".format(_figinfo["_fail_cutoffs"]["_riboScatter_cutoff"])  + ")",
                    "Warn (" + "{:.0%}".format(_figinfo["_warn_cutoffs"]["_riboScatter_cutoff"]) + ")",
                    "Mean rRNA/Aligned Reads (" + "{:.0%}".format(_slope_current) + ")"],
                loc='upper left',
                frameon=False,
                fontsize=_figinfo["_legend_size"])

    _ax = mk_axes(_ax) 
    _ax = needs_fail_or_warn(_ax,_slope_current,_figinfo,"_riboScatter_cutoff","upper")
   
    return _f

#### Plot 6: Sequence Contamination - Violin Plot ####
def plotViolin_dualAxis(_input_tup, _userDf, _background_df, _position,_figinfo,_f=None):
    # Load the sequence contamination levels for the background data
    _contaminant_df_untrim = _background_df[["Percent_Overrepresented_Seq_Untrimmed", "Percent_Adapter_Content_Untrimmed"]]
    _contaminant_df_trim = _background_df[["Percent_Overrepresented_Seq_Trimmed", "Percent_Adapter_Content_Trimmed"]]
    
    _contaminant_df_untrim.columns = ["Overrepresented", "Adapter"]
    _contaminant_df_trim.columns = ["Overrepresented", "Adapter"]

    _contaminant_melt_untrim = pd.melt(_contaminant_df_untrim, var_name="Contamination_Metric", value_name="Percent")
    _contaminant_melt_trim = pd.melt(_contaminant_df_trim, var_name="Contamination_Metric", value_name="Percent")

    _current_overrep_untrim = _input_tup[8] # if these are hard referencing column indices this has the potational to be a huge issue
    _current_adapter_untrim = _input_tup[9]
    _current_overrep_trim = _input_tup[10]
    _current_adapter_trim = _input_tup[11]
    
    # find the maximum value across each plot for each plot
    _max_array_untrim = _contaminant_df_untrim.max().values
    _max_array_untrim = np.append(_max_array_untrim,np.max(_current_overrep_untrim))
    _max_array_untrim = np.append(_max_array_untrim,np.max(_current_adapter_untrim))
    
    _max_array_trim = _contaminant_df_untrim.max().values
    _max_array_trim = np.append(_max_array_trim,np.max(_current_overrep_trim))
    _max_array_trim = np.append(_max_array_trim,np.max(_current_adapter_trim))
    
    # Remove the current batch mean from the USER dataframe
    _user_minusBatchMean_df = _userDf.drop(_userDf.tail(1).index)

    _mean_overrep_untrim = _user_minusBatchMean_df.loc[:, 'Percent_Overrepresented_Seq_Untrimmed'].mean()
    _mean_adapter_untrim = _user_minusBatchMean_df.loc[:, 'Percent_Adapter_Content_Untrimmed'].mean()

    _mean_overrep_trim = _user_minusBatchMean_df.loc[:, 'Percent_Overrepresented_Seq_Trimmed'].mean()
    _mean_adapter_trim = _user_minusBatchMean_df.loc[:, 'Percent_Adapter_Content_Trimmed'].mean()



    # Define color palette
    _contaminant_pal = {"Overrepresented": "lightgray", "Adapter": "gray"}

    _gridsp = matplotlib.gridspec.GridSpec(_figinfo["_subplot_rows"]*2, 2, figure=_f)
    
    _axis = _f.add_subplot(_gridsp[4, 1:])
    _axis2 = _f.add_subplot(_gridsp[5, 1:])

    sns.violinplot(x="Percent", y="Contamination_Metric", data=_contaminant_melt_untrim, palette=_contaminant_pal,
                   inner=None,ax=_axis, linewidth=0.3, orient="h", scale="count")

    sns.violinplot(x="Percent", y="Contamination_Metric", data=_contaminant_melt_trim, palette=_contaminant_pal,
                   inner=None,ax=_axis2, linewidth=0.3, orient="h", scale="count")

    _axis.set_title("Sequence Contamination", fontsize=_figinfo["_title_size"], pad=0)
    _axis  = set_ticks(_axis, _figinfo["_tick_size"])
    _axis2 = set_ticks(_axis2,_figinfo["_tick_size"])

    _axis.set_xlabel("\n\n\n\n\n\n\n  ", fontsize=_figinfo["_label_size"], labelpad=0.5)
    _axis.set_ylabel("")

    _axis2.set_xlabel("(%)", fontsize=_figinfo["_label_size"], labelpad=0.5)
    _axis2.set_ylabel("")

    _axis.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    _axis2.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())

    _axis.set_yticklabels(['Untrimmed \nOverrepresented', 'Untrimmed \nAdapter'])
    _axis2.set_yticklabels(['Trimmed \nOverrepresented', 'Trimmed \nAdapter'])

    
    _axis  = mk_axes(_axis)
    _axis2 = mk_axes(_axis2)

    _x_bottom, _x_top = _axis.get_xlim()
    _y_bottom, _y_top = _axis.get_ylim()

    ### Adding cutoff markers
    _axis3 = _axis2.twinx()

    _axis2,_axis3 = mk_axes(_axis2,_axis3)
    _axis3.yaxis.set_ticks([])

    _axis3.xaxis.label.set_visible(False)

    _axis3.set_ylim(_axis2.get_ylim()[0], _axis2.get_ylim()[1])

    _markers = [ # overrepresented
    (_axis3, _figinfo["_fail_cutoffs"]["_violin_cutoff_overrep_trimmed"], -0.4, 'Fail', _figinfo["_fail_color"]),
    (_axis3, _figinfo["_warn_cutoffs"]["_violin_cutoff_overrep_trimmed"], -0.6,'Warn', _figinfo["_warn_color"]),
    # adapter 
    (_axis3, _figinfo["_fail_cutoffs"]["_violin_cutoff_adapter_trimmed"], .75, 'Fail', _figinfo["_fail_color"]),
    (_axis3, _figinfo["_warn_cutoffs"]["_violin_cutoff_adapter_trimmed"], .85, 'Warn', _figinfo["_warn_color"])]

    for _axs, _cutoff, yloc,label, color in _markers:
        _axs.plot(_cutoff, yloc, marker='v', ms=1, c=color, clip_on=False)
        _axs.text(_cutoff, yloc - .1 , label, fontsize=_figinfo["_tick_size"], color=color, horizontalalignment='center')

    _line_overrep_untrim = _axis.axvline(x=_current_overrep_untrim, ymin=0.5, ymax=0.95, alpha=0.8, color=_figinfo["_curr_sample_color"],
                                         linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_overrep_untrim))
    _line_mean_overrep_untrim = _axis.axvline(x=_mean_overrep_untrim, ymin=0.5, ymax=0.95, alpha=0.8, color='indigo',
                                              linestyle='--', linewidth=0.35,
                                              label='{:.2f}%'.format(_mean_overrep_untrim))

    _line_adapter_untrim = _axis.axvline(x=_current_adapter_untrim, ymin=0.05, ymax=0.45, alpha=0.8, color=_figinfo["_curr_sample_color"],
                                         linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_adapter_untrim))

    _line_mean_adapter_untrim = _axis.axvline(x=_mean_adapter_untrim,ymin=.05,ymax=.45,alpha=.8, color ="indigo",
                                              linestyle='--', linewidth=0.35,
                                              label='{:.2f}%'.format(_mean_adapter_untrim))

    _line_overrep_trim = _axis2.axvline(x=_current_overrep_trim, ymin=0.55, ymax=0.95, alpha=0.8, color=_figinfo["_curr_sample_color"],
                                        linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_overrep_trim))
    _line_mean_overrep_trim = _axis2.axvline(x=_mean_overrep_trim, ymin=0.55, ymax=0.95, alpha=0.8, color='indigo',
                                             linestyle='--', linewidth=0.35, label='{:.2f}%'.format(_mean_overrep_trim))

    _line_adapter_trim = _axis2.axvline(x=_current_adapter_trim, ymin=0.05, ymax=0.45, alpha=0.8, color=_figinfo["_curr_sample_color"],
                                        linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_adapter_trim))
    _line_mean_adapter_trim = _axis2.axvline(x=_mean_adapter_trim, ymin=0.05, ymax=0.45, alpha=0.8, color='indigo',
                                             linestyle='--', linewidth=0.35, label='{:.2f}%'.format(_mean_adapter_trim))

    _axis.le

gend([_line_overrep_trim, _line_mean_overrep_trim], ["Current Sample", "Batch Mean"], loc='upper right',
                 frameon=False, ncol=1, fontsize=_figinfo["_legend_size"])

    plt.subplots_adjust(hspace=0)

    needs_fail_or_warn(_axis,_current_overrep_trim, _figinfo,"_violin_cutoff_overrep_trimmed","higher")
    needs_fail_or_warn(_axis,_current_overrep_trim, _figinfo,"_violin_cutoff_adapter_trimmed","higher") 
  #  if(_current_overrep_trim > _figinfo["_fail_violin_cutoff_overrep_trimmed"] or _current_adapter_trim > _figinfo["_fail_violin_cutoff_adapter_trimmed"]):
  #      insert_flag_fail(_axis)
  #  elif(_current_overrep_trim > _figinfo["_warn_violin_cutoff_overrep_trimmed"] or _current_adapter_trim > _figinfo["_warn_violin_cutoff_adapter_trimmed"]):
  #      insert_flag_warn(_axis)

    return _f


# function designed to return pvalues for GCinformation if supplied
def GC_KSstats(_coverage_df):
    
    # initialize list to store values
    _kslst = []
    # calculate mean GBC for the whole library
    _mean_df = pd.DataFrame()
    _mean_df["gc_mean"] = _coverage_df.median(axis=1)
    
    for column_name, _column_data in _coverage_df.iteritems():
        _ks_stat, _ks_pval = stats.ks_2samp(_column_data, _mean_df['gc_mean'])
        _kslst.append(_ks_stat)

    return _kslst

#  GeneBody Coverage Plot
def plotGC(_ipTuple, _coverage_df, _position, _plot_title,_figinfo,_fig=None):
    
    _axis = _fig.add_subplot(_figinfo["_subplot_rows"], 2, _position)

    # Calculate mean GeneBody Coverage for the entire library
    _mean_df = pd.DataFrame()
    _mean_df['gc_mean'] = _coverage_df.median(axis=1)

    # acquire the pvalue information from the _figinfo object
    _ks_pvals = _figinfo['_gbc_pvals']
    _ks_pval  = _ks_pvals[_ipTuple[0]]
    # Calculate Confidence Interval for the mean GC line

    # Calculate 95% interval for each position
    _err = _coverage_df.std(axis=1)*2

    # Plot current sample with library mean
    _x = np.arange(1, 101, 1)


    _axis.plot(_x, _coverage_df, color="lightgray",alpha = .4, linewidth=0.5, linestyle='-')
    _axis.plot(_x, _coverage_df[_ipTuple[1]], color=_figinfo["_curr_sample_color"], linewidth=0.5, linestyle='-')
    _axis.plot(_x, _mean_df['gc_mean'], color='indigo', linewidth=0.5, linestyle='--', alpha=0.8)

    _axis.fill_between(_x, _mean_df['gc_mean'] - _err, _mean_df['gc_mean'] + _err, facecolor='yellow', alpha=0.5)

    _axis = set_ticks(_axis,_figinfo["_tick_size"])

    _axis.set_xlim(0, 105)
    _axis.set_title(_plot_title, fontsize=_figinfo["_title_size"])

    _axis.set_xlabel("Gene Percentile (5' " + u"\u2192" + " 3')", fontsize=_figinfo["_label_size"], labelpad=2)
    _axis.set_ylabel("Read Density", fontsize=_figinfo["_label_size"], labelpad=2)
    
    # make the symbols for the legend
    _current_sample_line = matplotlib.lines.Line2D([0], [0], color=_figinfo["_curr_sample_color"], linewidth=0.5, linestyle='-', alpha=0.8)
    _background_lines = matplotlib.lines.Line2D([0], [0], color="lightgray", linewidth=0.5, linestyle='-', alpha=0.6)
    _library_line = matplotlib.lines.Line2D([0], [0], color="indigo", linewidth=0.5, linestyle='--', alpha=0.8)
    _extra_confidenceInterval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='yellow', fill=True,
                                                             edgecolor='yellow', linewidth=1.2, alpha=0.5)
    _extra_ksPval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None', linewidth=0)

    _axis.legend([_current_sample_line, _library_line,_background_lines, _extra_ksPval],
                 ["Current Sample", "Batch Mean",
                 "Batch Samples",  "KS Pvalue: " + str(round(_ks_pval, 3))],
                 loc='lower center', frameon=False, fontsize=_figinfo["_legend_size"], ncol=1)

    _axis = mk_axes(_axis)
    _axis = needs_fail_or_warn(_axis,_ks_pval,_figinfo,"_alpha","lower")

    return _fig

def calcHistPval(_hist_df):

    _data_df = _hist_df.drop(['Unnamed: 0'], axis=1)
    # code for calculating Z value of number of expressed genes. In need of some improvement.
    _sum_df = _data_df.sum().round()
    _zscore = stats.zscore(_sum_df)
    _pvals  = stats.norm.sf(abs(_zscore))
    
    return(_pvals)

# Plot 8 : Gene Expression Distribution Plot 
def plotNegBin(_ipTuple, _hist_df, _user_df,_pos, _plot_title,_figinfo,_f=None):
    _index_array = _hist_df.iloc[:, 0]

    _low_vals = []
    _high_vals = []

    for _i in _index_array:
        _low_vals.append(float(_i.strip('(').strip(']').split(',')[0]))
        _high_vals.append(float(_i.strip('(').strip(']').split(',')[1]))

    _x_vals = _low_vals

    ## Preparing the data_df and libMean_df for all bins
    _data_df = _hist_df.drop(['Unnamed: 0'], axis=1)
    _libMean_df = pd.DataFrame()
    _libMean_df['Mean'] = _data_df.iloc[:, :-1].mean(numeric_only=True, axis=1)

    _data_df_dropped = _data_df
    _mean_df_dropped = _libMean_df

    _max_df = pd.DataFrame()
    _max_df['max_val'] = _data_df.max(axis=1)
    _max_df.index == _ipTuple[1]
    _mean_array = _mean_df_dropped.Mean.values
    _current_samp_array = _data_df_dropped[_ipTuple[1]].values

    # code for calculating Z value of number of expressed genes. In need of some improvement.
    _sum_df = _data_df_dropped.sum().round()
    _curr_sum = _current_samp_array.sum().round()
    _curr_ndx = np.where(_sum_df == _curr_sum)[0][0]
    _zscore = stats.zscore(_sum_df)
    _pvals  = stats.norm.sf(abs(_zscore))
    _curr_pval = _pvals[_curr_ndx]

    _ax = _f.add_subplot(_figinfo["_subplot_rows"], 2, _pos)

    _col_names = [_cl for _cl in _data_df_dropped.columns]

    ## Plotting all current library distributions with the current sample highlighted on each page
    for _col in _col_names:

        if _col == _ipTuple[1]:
            plt.plot(_x_vals, _data_df_dropped[_col], color=_figinfo["_curr_sample_color"], linewidth=0.5, linestyle='-', zorder=24)
        else:
            plt.plot(_x_vals, _data_df_dropped[_col], color='silver', linewidth=0.5, linestyle='-')

    ## Plotting the mean distribution
    _ax.plot(_x_vals, _mean_df_dropped['Mean'], color='indigo', linewidth=0.5, linestyle='--', alpha=0.8, zorder=23)

    _ax = set_ticks(_ax,_figinfo["_tick_size"])

    _ax.set_xlim(0, 10)
    _ax.set_ylim(0, _max_df['max_val'].max())

    _ax.set_title(_plot_title, fontsize=_figinfo["_title_size"])

    _ax.set_xlabel("Expression Level (log2(CPM)+1)", fontsize=_figinfo["_label_size"], labelpad=2)
    _ax.set_ylabel("Frequency", fontsize=_figinfo["_label_size"], labelpad=2)

    _ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(_x_vals[::5]))
    _ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(_x_vals[::5]))
    _ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(500))
    _ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))

    _current_samp_line = matplotlib.lines.Line2D([0], [0], color=_figinfo["_curr_sample_color"], linewidth=0.5, linestyle='-', alpha=0.8)
    _lib_line = matplotlib.lines.Line2D([0], [0], color="indigo", linewidth=0.5, linestyle='--', alpha=0.8)

    _extra_Ztest_stat = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None',
                                                     linewidth=0)
    _extra_Ztest_Pval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None',
                                                     linewidth=0)

    _ax.legend([_current_samp_line, _lib_line, _extra_Ztest_Pval],
               ["Current Sample", "Batch  Mean", "Pvalue (# Detected Genes): " + str(round(_curr_pval.item(), 3))], loc='upper right',
               frameon=False, fontsize=_figinfo["_legend_size"], ncol=1)
    _ax = mk_axes(_ax)
    _ax = needs_fail_or_warn(_ax,_curr_pval,_figinfo,"_alpha","lower")

    return _f

if __name__ == "__main__":
    main()
