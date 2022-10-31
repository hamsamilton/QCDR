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
#from pandas.plotting import register_matplotlib_converters
from statsmodels.stats.weightstats import ztest



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





####### Add all the data loading, cleaning and other helper functions here ##########

def fmt_number(number, pos=None):
    if number == 0:
        return '0M'

    else:

        magnitude = 0
        while abs(number) >= 1000:
            magnitude += 1
            number /= 1000.0
        return '%.0f%s' % (number, ['', 'K', 'M', 'B', 'T', 'Q'][magnitude])


def fmt_scatter_million(_x, _pos):
    return '%1.1f' % (_x * 1e-6)


def fmt_interval(_x, _pos):
    return '0:.0s%'.format(_x)


def fmt_million(_x, _pos):
    # return '{0:.0f}M'.format(_x)
    return '{0:.0f}'.format(_x)


def fmt_contaminant(_c, _pos):
    return '{0:.2f}%'.format(_c)


def fmt(_x, _pos):
    # return '{0:.0f}%'.format(_x)
    return '{0:.0f}'.format(_x)


def fmt_cov(_x, _pos):
    return '{0:.2f}'.format(_x)


def byMillion(_lib):
    return _lib / 1000000


def insert_flag_image(_image, loc=3, ax=None, zoom=1, **kw):
    if ax == None:
        ax = plt.gca()

    _im_box = matplotlib.offsetbox.OffsetImage(_image, zoom=zoom * 0.05)
    _anch_box = matplotlib.offsetbox.AnchoredOffsetbox(loc=loc, child=_im_box, frameon=False, **kw)

    ax.add_artist(_anch_box)


def insert_flag_fail(ax=None):
    if ax == None:
        ax = plt.gca()

    anch_text = matplotlib.offsetbox.AnchoredText("FAILURE ", pad=0.0001, borderpad=1,
                                                  loc=3,
                                                  prop=dict(size=4, snap=True, backgroundcolor='tomato', color='yellow',
                                                            alpha=0.9, zorder=5, fontweight='roman', fontfamily='serif',
                                                            fontstretch='extra-expanded'),
                                                  frameon=True,
                                                  bbox_to_anchor=(0., 1.03),
                                                  bbox_transform=ax.transAxes)

    # anch_text.patch.set_boxstyle=("round")

    ax.add_artist(anch_text)

    return None


def insert_flag_warn(ax=None):
    if ax == None:
        ax = plt.gca()

    anch_text = matplotlib.offsetbox.AnchoredText("WARNING", pad=0.00000001, borderpad=1,
                                                  loc=3, prop=dict(size=3.7, snap=True, animated=True,
                                                                   backgroundcolor='yellow', color='red', alpha=0.9,
                                                                   zorder=5, fontweight='roman', fontfamily='serif',
                                                                   fontstretch='expanded'),
                                                  frameon=True,
                                                  bbox_to_anchor=(0., 1.03),
                                                  bbox_transform=ax.transAxes)

    # anch_text.patch.set_boxstyle=("round")

    ax.add_artist(anch_text)

    return None


def ztest_prob(dist_current, dist_mean, val):
    # print(dist_current)
    # print(dist_mean)

    ztest_stat, ztest_pval = ztest(dist_current, dist_mean, value=val, alternative='two-sided')

    # print(ztest_stat, ztest_pval)

    return ztest_stat, ztest_pval


def cdf_prob(df):


    df_cdf = df.apply(lambda x: 1 - stats.norm.cdf(x))
    df_pvalue = df.apply(lambda x: stats.norm.sf(x))
    print(df_pvalue)

    return df_cdf




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
    # xytext = (10, 10)
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




####### Add all the plotting functions here for all panels ########


def plotHist_ipSize_old(_in_tuple, _slaveDf, _master_df, _pos, _f=None):
    # print(_master_df.Input_Size)
    bin = np.arange(0, 100 + 1, 1)

    print(_master_df)
    print(_master_df.columns)

    ### Add a hard-clip limit of 40 million reads to the master['Input_Size']
    _master_df.loc[:, 'Input_Size'] = _master_df.loc[:, 'Input_Size'].clip(upper=40)
    print(_master_df.Input_Size)

    print(_in_tuple)

    _slaveDf.loc[:, "Input_Size"] = _slaveDf.loc[:, "Input_Size"].clip(upper=40000000)
    print(_slaveDf)
    print(_slaveDf.Input_Size)

    _out, _bins = pd.cut(_master_df['Input_Size'], bins=bin, retbins=True, right=True, include_lowest=False)
    print(_out)
    print(_bins)
    _xlabs = [str(xt) for xt in _bins[0::5]]

    print(_in_tuple)

    _ip_norm = _slaveDf.loc[:, 'Input_Size'].apply(byMillion)
    # print(_ip_norm[_in_tuple.Index])

    _lib_mean = _ip_norm.mean()
    _current_sample = _ip_norm[_in_tuple.Index]

    # print(_lib_mean)

    if not _f is None:
        plt.gcf()

    _ax = _f.add_subplot(3, 2, _pos)

    # _ax = _master_df['Input_Size'].plot(kind='hist', bins=_bins, ax=plt.gca())
    # _ax1 = _master_df.plot(y='Input_Size', kind='kde', secondary_y=True, mark_right=True, legend=False, lw=0.5, ax=_ax, color='magenta')

    # _ax.hist(_master_df['Input_Size'], bins=_bins)
    _master_df["Input_Size"].plot(kind='hist', bins=_bins, ax=_ax, color='lightgray')

    _ax1 = _ax.twinx()
    _master_df.plot(y="Input_Size", kind='kde', legend=False, ax=_ax1, color='dimgray', mark_right=True, lw=0.7,
                    alpha=0.8)

    _ax.tick_params(axis='x', which='both', length=1, width=0.5, labelbottom=True, bottom=True, labelsize=3,
                    direction='out', pad=2)
    _ax.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=4, labelleft=True, left=True,
                    direction='out', pad=2)

    _ax.set_xlim(0, 40)

    _ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(_bins[0::5]))
    _ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_million))

    _ax.set_title("Sequencing Depth", fontsize=4.4)
    _ax.set_xlabel('Total Reads (Millions)', labelpad=1, fontsize=3)

    _ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    _ax.set_ylabel('Frequency', labelpad=2, fontsize=3)

    _ax1.yaxis.set_ticks([])
    _ax1.yaxis.label.set_visible(False)
    # _ax1.set_ylabel("Kernel Density Estimate", fontsize=4, labelpad=3)

    _ax.spines['top'].set_visible(False)
    _ax.spines['right'].set_visible(False)
    _ax.spines['left'].set_visible(True)
    _ax.spines['bottom'].set_visible(True)
    _ax.spines['left'].set_color('black')
    _ax.spines['bottom'].set_color('black')

    _ax.spines['left'].set_linewidth(0.55)
    _ax.spines['bottom'].set_linewidth(0.55)

    # _ax1.spines['top'].set_visible(False)
    _ax1.spines['right'].set_visible(False)
    _ax1.spines['bottom'].set_visible(False)
    _ax1.spines['left'].set_visible(False)

    if _current_sample > _lib_mean:
        _adj_flag = True
        _rotation_current = 270
        _rotation_mean = 90
    elif _current_sample < _lib_mean:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270
    else:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270

    ### Adding cutoff markers
    _ax.plot(_ipReads_cutoff_fail, _ax.get_ylim()[1] - 1, marker='v', ms=0.8, c='red')
    _ax.text(_ipReads_cutoff_fail, _ax.get_ylim()[1] - 0.4, 'Fail', fontsize=2, color='red',
             horizontalalignment='center')

    _ax.plot(_ipReads_cutoff_warn, _ax.get_ylim()[1] - 1, marker='v', ms=0.8, c='gold')
    _ax.text(_ipReads_cutoff_warn, _ax.get_ylim()[1] - 0.4, 'Warn', fontsize=2, color='gold',
             horizontalalignment='center')

    # Current Sample Line and Label
    _line1 = _ax.axvline(x=_current_sample, alpha=0.8, color='red', linestyle='-', linewidth=0.5,
                         label='{:2.2f}M'.format(_current_sample))

    if _adj_flag:
        _ax.text(_current_sample + 0.5, (_ax.get_ylim()[1] / 2), '{:2.2f}M'.format(_current_sample),
                 rotation=_rotation_current, fontsize=2.5, zorder=2)
    else:
        _ax.text(_current_sample - 2, (_ax.get_ylim()[1] / 2), '{:2.2f}M'.format(_current_sample),
                 rotation=_rotation_current, fontsize=2.5, zorder=2)

    # Current Library Mean Line and Label
    _line2 = _ax.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5,
                         label='{:2.2f}M'.format(_lib_mean))
    if _adj_flag:
        _ax.text(_lib_mean - 2, ((_ax.get_ylim()[1] / 2) + 1), '{:2.2f}M'.format(_lib_mean), rotation=_rotation_mean,
                 fontsize=2.5, zorder=2)
    else:
        _ax.text(_lib_mean + 0.5, ((_ax.get_ylim()[1] / 2) + 1), '{:2.2f}M'.format(_lib_mean), rotation=_rotation_mean,
                 fontsize=2.5, zorder=2)

    _kde_line = matplotlib.lines.Line2D([0], [0], color="gray", linewidth=0.5, linestyle='-')
    _ax.legend([_line1, _line2], ["Current Sample", "Library Mean"], loc='best', frameon=False, fontsize=3)

    _ax.set_facecolor('white')

    # Flagging for FAILURE/WARNING
    if _current_sample <= _ipReads_cutoff_fail:
        print("Sample flagged for FAILURE!!!")
        insert_flag_fail(_ax)

    elif (_current_sample <= _ipReads_cutoff_warn) and (_current_sample > _ipReads_cutoff_fail):
        print("Sample flagged for WARNING!!!")
        insert_flag_warn(_ax)

    else:
        print("Sample PASSED cutoffs!")

    plt.box(on=None)

    return _f


#### Plot 1: Input Size ####
def plotHist_ipSize(_in_tuple, _userDf, _background_df, _pos, _f=None):
    # print(_background_df.Input_Size)
    bin = np.arange(0, 100 + 1, 1)

    print(_background_df)
    print(_background_df.columns)

    ### Add a hard-clip limit of 60 million reads to the master['Input_Size']
    _background_df.loc[:, 'Input_Size'] = _background_df.loc[:, 'Input_Size'].clip(upper=40)
    print(_background_df.Input_Size)

    ### CDF cutoff
    #cdf_prob(_background_df.Input_Size)




    print(_in_tuple)

    _userDf.loc[:, "Input_Size"] = _userDf.loc[:, "Input_Size"].clip(upper=40000000)
    print(_userDf)
    print(_userDf.Input_Size)

    _out, _bins = pd.cut(_background_df['Input_Size'], bins=bin, retbins=True, right=True, include_lowest=False)
    print(_out)
    print(_bins)
    _xlabs = [str(xt) for xt in _bins[0::5]]

    print(_in_tuple)

    _ip_norm = _userDf.loc[:, 'Input_Size'].apply(byMillion)
    # print(_ip_norm[_in_tuple.Index])

    _lib_mean = _ip_norm.mean()
    _current_sample = _ip_norm[_in_tuple.Index]

    # print(_lib_mean)

    if not _f is None:
        plt.gcf()

    _ax = _f.add_subplot(3, 2, _pos)

    # _ax = _background_df['Input_Size'].plot(kind='hist', bins=_bins, ax=plt.gca())
    # _ax1 = _background_df.plot(y='Input_Size', kind='kde', secondary_y=True, mark_right=True, legend=False, lw=0.5, ax=_ax, color='magenta')

    # _ax.hist(_background_df['Input_Size'], bins=_bins)
    _background_df["Input_Size"].plot(kind='hist', bins=_bins, ax=_ax, color='lightgray')

    _ax1 = _ax.twinx()
    #_background_df.plot(y="Input_Size", kind='kde', legend=False, ax=_ax1, color='dimgray', mark_right=True, lw=0.7, alpha=0.8)
    sns.distplot(_background_df["Input_Size"], hist=False, bins=_bins, ax=_ax1, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})

    _ax.tick_params(axis='x', which='both', length=1, width=0.5, labelbottom=True, bottom=True, labelsize=3, direction='out', pad=2)
    _ax.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=4, labelleft=True, left=True, direction='out', pad=2)

    _ax.set_xlim(0, 40)

    _ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(_bins[0::5]))
    _ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_million))

    _ax.set_title("Sequencing Depth", fontsize=4.4)
    _ax.set_xlabel('Total Reads (Millions)', labelpad=1, fontsize=3.5)

    _ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    _ax.set_ylabel('Frequency', labelpad=2, fontsize=3.5)

    _ax.set_facecolor('white')

    _ax1.yaxis.set_ticks([])
    _ax1.yaxis.label.set_visible(False)
    # _ax1.set_ylabel("Kernel Density Estimate", fontsize=4, labelpad=3)

    _ax.spines['top'].set_visible(False)
    _ax.spines['right'].set_visible(False)
    _ax.spines['left'].set_visible(True)
    _ax.spines['bottom'].set_visible(True)
    _ax.spines['left'].set_color('black')
    _ax.spines['bottom'].set_color('black')

    _ax.spines['left'].set_linewidth(0.55)
    _ax.spines['bottom'].set_linewidth(0.55)

    # _ax1.spines['top'].set_visible(False)
    _ax1.spines['right'].set_visible(False)
    _ax1.spines['bottom'].set_visible(False)
    _ax1.spines['left'].set_visible(False)

    if _current_sample > _lib_mean:
        _adj_flag = True
        _rotation_current = 270
        _rotation_mean = 90
    elif _current_sample < _lib_mean:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270
    else:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270

    ### Adding cutoff markers
    _ax.plot(_ipReads_cutoff_fail, _ax.get_ylim()[1] - 1, marker='v', ms=0.8, c='red')
    _ax.text(_ipReads_cutoff_fail, _ax.get_ylim()[1] - 0.4, 'Fail', fontsize=2, color='red',
             horizontalalignment='center')

    _ax.plot(_ipReads_cutoff_warn, _ax.get_ylim()[1] - 1, marker='v', ms=0.8, c='gold')
    _ax.text(_ipReads_cutoff_warn, _ax.get_ylim()[1] - 0.4, 'Warn', fontsize=2, color='gold',
             horizontalalignment='center')

    # Current Sample Line and Label
    _line1 = _ax.axvline(x=_current_sample, alpha=0.8, color='red', linestyle='-', linewidth=0.5,
                         label='{:2.2f}M'.format(_current_sample))

    if _adj_flag:
        _ax.text(_current_sample + 0.5, (_ax.get_ylim()[1] / 2), '{:2.2f}M'.format(_current_sample),
                 rotation=_rotation_current, fontsize=2.5, zorder=2)
    else:
        _ax.text(_current_sample - 2, (_ax.get_ylim()[1] / 2), '{:2.2f}M'.format(_current_sample),
                 rotation=_rotation_current, fontsize=2.5, zorder=2)

    # Current Library Mean Line and Label
    _line2 = _ax.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5,
                         label='{:2.2f}M'.format(_lib_mean))
    if _adj_flag:
        _ax.text(_lib_mean - 2, ((_ax.get_ylim()[1] / 2) + 1), '{:2.2f}M'.format(_lib_mean), rotation=_rotation_mean,
                 fontsize=2.5, zorder=2)
    else:
        _ax.text(_lib_mean + 0.5, ((_ax.get_ylim()[1] / 2) + 1), '{:2.2f}M'.format(_lib_mean), rotation=_rotation_mean,
                 fontsize=2.5, zorder=2)

    _kde_line = matplotlib.lines.Line2D([0], [0], color="gray", linewidth=0.5, linestyle='-')
    _ax.legend([_line1, _line2], ["Current Sample", "Batch Mean"], loc='best', frameon=False, fontsize=3)



    # Flagging for FAILURE/WARNING
    if _current_sample <= _ipReads_cutoff_fail:
        print("Sample flagged for FAILURE!!!")
        insert_flag_fail(_ax)

    elif (_current_sample <= _ipReads_cutoff_warn) and (_current_sample > _ipReads_cutoff_fail):
        print("Sample flagged for WARNING!!!")
        insert_flag_warn(_ax)

    else:
        print("Sample PASSED cutoffs!")

    plt.box(on=None)


    return _f



#### Plot 2 : Trimming Percentage ####
def plotHist_trimming(_ip_tuple, _user_df, _retro_df, _colname, _plot_label, _position, _figure=None):
    
    
    bin_data = np.arange(0, 100 + 1, 1)

    _retro_df.loc[:, "Percent_PostTrim"] = _retro_df.loc[:, "Percent_PostTrim"].clip(lower=60)
    _user_df.loc[:, "Percent_PostTrim"] = _user_df.loc[:, "Percent_PostTrim"].clip(lower=60)

    # print(_retro_df.Percent_PostTrim)
    # print(_ip_tuple[3])
    #print(_user_df.Percent_PostTrim)

    _out, _bins = pd.cut(_retro_df[_colname], bins=bin_data, retbins=True, right=True, include_lowest=True)
    # print(_out.value_counts(sort=True).index.categories.tolist())
    # print(_bins)
    _xtick_labels = pd.Series(_out.value_counts(sort=True).index.categories)
    # LabelEncoder().fit([x for y in _xtick_labels.get_values() for x in y])
    # print(_xtick_labels.get_values())

    _xtick_labs = [str(xt) for xt in _bins[0::5]]
    # _xtick_labs = _xtick_labs.tolist()

    _user_minusBatchMean_df = _user_df.drop(_user_df.tail(1).index)

    #print(_user_minusBatchMean_df)
    
    

    if _colname == "Percent_PostTrim":

        _current_sample = _ip_tuple.Percent_PostTrim

        if _current_sample < 60:
            _current_sample = 60
        else:
            pass

        _lib_mean = _user_minusBatchMean_df[_colname].mean()

    else:
        _current_sample = 0
        _lib_mean = 0
        print("Not a legitimate metric!\n")


    #print(_ip_tuple)
    #print(_current_sample)
    #print(_lib_mean)



    if not _figure is None:
        plt.gcf()

    _axis = _figure.add_subplot(3, 2, _position)

    # _axis = _retro_df[_colname].plot(kind='hist', bins=_bins, ax=_axis)
    # _axis1 = _retro_df.plot(y=_colname, kind='kde', secondary_y=True, mark_right=False, legend=False, lw=0.4,color='magenta', ax=_axis)

    _axis.hist(x=_retro_df[_colname], bins=_bins, histtype='bar', color='lightgray')

    _axis1 = _axis.twinx()
    #_retro_df.plot(y=_colname, kind = 'kde', legend = False, ax = _axis1, color = 'dimgray', mark_right = True, lw = 0.7, alpha = 0.8)
    sns.distplot(_retro_df["Percent_PostTrim"], hist=False, bins=_bins, ax=_axis1, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})


    _axis.set_xlim(59, 103)

    _axis.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(_bins[0::5]))
    _axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))

    _axis.tick_params(axis='x', which='both', length=1, width=0.5, labelbottom=True, bottom=True, labelsize=3,
                      direction='out', pad=2)
    _axis.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=4, labelleft=True, left=True,
                      direction='out', pad=2)

    _axis.set_title(_plot_label, fontsize=4.4)

    if _colname == "Percent_PostTrim":
        _axis.set_xlabel('% Post-Trim / Total Reads', labelpad=1, fontsize=3)

        if _current_sample <= _trimmedReads_cutoff_fail:
            print("Sample flagged for FAILURE!!!")
            insert_flag_fail(_axis)

        elif (_current_sample > _trimmedReads_cutoff_fail) and (_current_sample <= _trimmedReads_cutoff_warn):
            print("Sample flagged for WARNING!!!")
            insert_flag_warn(_axis)

        else:
            print("Sample PASSED cutoffs!")
            pass

        ### Adding cutoff markers
        _axis.plot(_trimmedReads_cutoff_fail, _axis.get_ylim()[1] - (_axis.get_ylim()[1] / 10), marker='v', ms=0.8,
                   c='red')
        _axis.text(_trimmedReads_cutoff_fail, _axis.get_ylim()[1] - (_axis.get_ylim()[1] / 20), 'Fail', fontsize=2,
                   color='red', horizontalalignment='center')

        _axis.plot(_trimmedReads_cutoff_warn, _axis.get_ylim()[1] - (_axis.get_ylim()[1] / 10), marker='v', ms=0.8,
                   c='gold')
        _axis.text(_trimmedReads_cutoff_warn, _axis.get_ylim()[1] - (_axis.get_ylim()[1] / 20), 'Warn', fontsize=2,
                   color='gold', horizontalalignment='center')




    
    _axis.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    _axis.set_ylabel('Frequency', labelpad=2, fontsize=3)

    # SECONDARY AXIS FOR KDE
    _axis1.get_yaxis().set_ticks([])
    _axis1.yaxis.label.set_visible(False)
    # _axis1.set_ylabel('Kernel Density Estimate', fontsize=4, labelpad=3)

    _axis.spines['top'].set_visible(False)
    _axis.spines['right'].set_visible(False)
    _axis.spines['left'].set_visible(True)
    _axis.spines['bottom'].set_visible(True)
    _axis.spines['left'].set_color('black')
    _axis.spines['bottom'].set_color('black')
    _axis.spines['bottom'].set_linewidth(0.55)
    _axis.spines['left'].set_linewidth(0.55)

    _axis1.spines['top'].set_visible(False)
    _axis1.spines['right'].set_visible(False)
    _axis1.spines['left'].set_visible(False)
    _axis1.spines['bottom'].set_visible(False)

    # _axis.grid(b=False)
    # _axis1.grid(b=False)

    if _current_sample > _lib_mean:
        _adj_flag = True
        _rotation_current = 270
        _rotation_mean = 90
    elif _current_sample < _lib_mean:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270
    else:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270

    # Current Sample Line and Label
    _line1 = _axis.axvline(x=_current_sample, alpha=0.8, color='red', linestyle='-', linewidth=0.5,
                           label='{:2.2f}%'.format(_current_sample))

    if _adj_flag:
        _axis.text(_current_sample + 0.5, (_axis.get_ylim()[1] / 2), '{:2.2f}%'.format(_current_sample),
                   rotation=_rotation_current, fontsize=2.5, zorder=2)
    else:
        _axis.text(_current_sample - 2, (_axis.get_ylim()[1] / 2), '{:2.2f}%'.format(_current_sample),
                   rotation=_rotation_current, fontsize=2.5, zorder=2)

    # Current Library Mean Line and Label
    _line2 = _axis.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5,
                           label='{:2.2f}%'.format(_lib_mean))

    if _adj_flag:
        _axis.text(_lib_mean - 2, ((_axis.get_ylim()[1] / 2) + 1), '{:2.2f}%'.format(_lib_mean),
                   rotation=_rotation_mean, fontsize=2.5, zorder=2)
    else:
        _axis.text(_lib_mean + 0.5, ((_axis.get_ylim()[1] / 2) + 1), '{:2.2f}%'.format(_lib_mean),
                   rotation=_rotation_mean, fontsize=2.5, zorder=2)

    ## Superimpose the Kernel Density Estimate line over the distribution
    _kde_line = matplotlib.lines.Line2D([0], [0], color="dimgray", linewidth=0.5, linestyle='-')

    _axis.legend([_line1, _line2], ["Current Sample", "Batch Mean"], loc='best', frameon=False, fontsize=3, ncol=1)
    _axis.set_facecolor('white')

    return _figure


#### Plot 3: Alignment Percentage ####
def plotHist_alignment(_ip_tuple, _user_df, _retro_df, _colname, _plot_label, _position, _figure=None):
    bin_data = np.arange(0, 100 + 1, 1)

    _retro_df.loc[:, "Percent_PostTrim"] = _retro_df.loc[:, "Percent_PostTrim"].clip(lower=60)
    _user_df.loc[:, "Percent_PostTrim"] = _user_df.loc[:, "Percent_PostTrim"].clip(lower=60)

    _out, _bins = pd.cut(_retro_df[_colname], bins=bin_data, retbins=True, right=True, include_lowest=True)
    # print(_out.value_counts(sort=True).index.categories.tolist())
    # print(_bins)
    _xtick_labels = pd.Series(_out.value_counts(sort=True).index.categories)
    # LabelEncoder().fit([x for y in _xtick_labels.get_values() for x in y])
    # print(_xtick_labels.get_values())

    _xtick_labs = [str(xt) for xt in _bins[0::5]]
    # _xtick_labs = _xtick_labs.tolist()

    _user_minusBatchMean_df = _user_df.drop(_user_df.tail(1).index)

    #print(_user_minusBatchMean_df)


    if _colname == "Percent_Uniquely_Aligned":
        _current_sample = _ip_tuple.Percent_Uniquely_Aligned
        _lib_mean = _user_minusBatchMean_df[_colname].mean()
        #print(_current_sample)
        #print(_lib_mean)

    else:
        _current_sample = 0
        _lib_mean = 0
        print("Not a legitimate metric!\n")

    if not _figure is None:
        plt.gcf()

    _axis_plt3 = _figure.add_subplot(3, 2, _position)

    # _axis_plt3 = _retro_df[_colname].plot(kind='hist', bins=_bins, ax=_axis_plt3)
    # _axis_plt31 = _retro_df.plot(y=_colname, kind='kde', secondary_y=True, mark_right=False, legend=False, lw=0.4,color='magenta', ax=_axis_plt3)

    _axis_plt3.hist(x=_retro_df[_colname], bins=_bins, histtype='bar', color='lightgray')

    _axis1_plt3 = _axis_plt3.twinx()
    #_retro_df.plot(y=_colname, kind='kde', legend=False, ax=_axis1_plt3, color='dimgray', mark_right=True, lw=0.7, alpha=0.8)
    sns.distplot(_retro_df[_colname], hist=False, bins=_bins, ax=_axis1_plt3, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})


    _axis_plt3.set_xlim(0, 105)

    _axis_plt3.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(_bins[0::5], nbins=21))
    _axis_plt3.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))

    _axis_plt3.tick_params(axis='x', which='both', length=1, width=0.5, labelbottom=True, bottom=True, labelsize=3, direction='out', pad=2)
    _axis_plt3.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=4, labelleft=True, left=True, direction='out', pad=2)

    _axis_plt3.set_title(_plot_label, fontsize=4.4)




    _axis_plt3.set_xlabel('% Uniquely Aligned / Post-Trim Reads', labelpad=1, fontsize=3)

    if _current_sample <= _uniqAligned_cutoff_fail:
        print("Sample flagged for FAILURE!!!")
        insert_flag_fail(_axis_plt3)

    elif (_current_sample > _uniqAligned_cutoff_fail) and (_current_sample <= _uniqAligned_cutoff_warn):
        print("Sample flagged for WARNING!!!")
        insert_flag_warn(_axis_plt3)

    else:
        print("Sample PASSED cutoffs!")
        pass

    ### Adding cutoff markers
    _axis_plt3.plot(_uniqAligned_cutoff_fail, _axis_plt3.get_ylim()[1] - (_axis_plt3.get_ylim()[1] / 10), marker='v', ms=0.8, c='red')
    _axis_plt3.text(_uniqAligned_cutoff_fail, _axis_plt3.get_ylim()[1] - (_axis_plt3.get_ylim()[1] / 20), 'Fail', fontsize=2, color='red', horizontalalignment='center')

    _axis_plt3.plot(_uniqAligned_cutoff_warn, _axis_plt3.get_ylim()[1] - (_axis_plt3.get_ylim()[1] / 10), marker='v', ms=0.8, c='gold')
    _axis_plt3.text(_uniqAligned_cutoff_warn, _axis_plt3.get_ylim()[1] - (_axis_plt3.get_ylim()[1] / 20), 'Warn', fontsize=2, color='gold', horizontalalignment='center')



    _axis_plt3.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    _axis_plt3.set_ylabel('Frequency', labelpad=2, fontsize=3)

    # SECONDARY AXIS FOR KDE
    _axis1_plt3.get_yaxis().set_ticks([])
    _axis1_plt3.yaxis.label.set_visible(False)
    # _axis1_plt3.set_ylabel('Kernel Density Estimate', fontsize=4, labelpad=3)

    _axis_plt3.spines['top'].set_visible(False)
    _axis_plt3.spines['right'].set_visible(False)
    _axis_plt3.spines['left'].set_visible(True)
    _axis_plt3.spines['bottom'].set_visible(True)
    _axis_plt3.spines['left'].set_color('black')
    _axis_plt3.spines['bottom'].set_color('black')
    _axis_plt3.spines['bottom'].set_linewidth(0.55)
    _axis_plt3.spines['left'].set_linewidth(0.55)

    _axis1_plt3.spines['top'].set_visible(False)
    _axis1_plt3.spines['right'].set_visible(False)
    _axis1_plt3.spines['left'].set_visible(False)
    _axis1_plt3.spines['bottom'].set_visible(False)

    # _axis_plt3.grid(b=False)
    # _axis1_plt3.grid(b=False)

    if _current_sample > _lib_mean:
        _adj_flag = True
        _rotation_current = 270
        _rotation_mean = 90
    elif _current_sample < _lib_mean:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270
    else:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270

    # Current Sample Line and Label
    _line1 = _axis_plt3.axvline(x=_current_sample, alpha=0.8, color='red', linestyle='-', linewidth=0.5, label='{:2.2f}%'.format(_current_sample))

    if _adj_flag:
        _axis_plt3.text(_current_sample + 0.5, (_axis_plt3.get_ylim()[1] / 2), '{:2.2f}%'.format(_current_sample), rotation=_rotation_current, fontsize=2.5, zorder=2)
    else:
        _axis_plt3.text(_current_sample - 2, (_axis_plt3.get_ylim()[1] / 2), '{:2.2f}%'.format(_current_sample), rotation=_rotation_current, fontsize=2.5, zorder=2)

    # Current Library Mean Line and Label
    _line2 = _axis_plt3.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5, label='{:2.2f}%'.format(_lib_mean))

    if _adj_flag:
        _axis_plt3.text(_lib_mean - 2, ((_axis_plt3.get_ylim()[1] / 2) + 1), '{:2.2f}%'.format(_lib_mean), rotation=_rotation_mean, fontsize=2.5, zorder=2)
    else:
        _axis_plt3.text(_lib_mean + 0.5, ((_axis_plt3.get_ylim()[1] / 2) + 1), '{:2.2f}%'.format(_lib_mean), rotation=_rotation_mean, fontsize=2.5, zorder=2)

    ## Superimpose the Kernel Density Estimate line over the distribution
    _kde_line = matplotlib.lines.Line2D([0], [0], color="dimgray", linewidth=0.5, linestyle='-')

    _axis_plt3.legend([_line1, _line2], ["Current Sample", "Batch Mean"], loc='best', frameon=False, fontsize=3, ncol=1)
    _axis_plt3.set_facecolor('white')

    return _figure



#### Plot 4: Gene Exon Mapping ####
def plotHist_exonMapping(_ip_tuple, _user_df, _retro_df, _colname, _plot_label, _position, _figure=None):
    bin_data = np.arange(0, 100 + 1, 1)

    _retro_df.loc[:, "Percent_PostTrim"] = _retro_df.loc[:, "Percent_PostTrim"].clip(lower=60)
    _user_df.loc[:, "Percent_PostTrim"] = _user_df.loc[:, "Percent_PostTrim"].clip(lower=60)

    _out, _bins = pd.cut(_retro_df[_colname], bins=bin_data, retbins=True, right=True, include_lowest=True)
    # print(_out.value_counts(sort=True).index.categories.tolist())
    # print(_bins)
    _xtick_labels = pd.Series(_out.value_counts(sort=True).index.categories)
    # LabelEncoder().fit([x for y in _xtick_labels.get_values() for x in y])
    # print(_xtick_labels.get_values())

    _xtick_labs = [str(xt) for xt in _bins[0::5]]
    # _xtick_labs = _xtick_labs.tolist()

    _user_minusBatchMean_df = _user_df.drop(_user_df.tail(1).index)

    #print(_user_minusBatchMean_df)

    if _colname == "Percent_Exonic":
        _current_sample = _ip_tuple.Percent_Exonic
        _lib_mean = _user_minusBatchMean_df[_colname].mean()
        #print(_current_sample)
        #print(_lib_mean)

    else:
        _current_sample = 0
        _lib_mean = 0
        print("Not a legitimate metric!\n")

    if not _figure is None:
        plt.gcf()

    _axis_plt4 = _figure.add_subplot(3, 2, _position)

    # _axis_plt4 = _retro_df[_colname].plot(kind='hist', bins=_bins, ax=_axis_plt4)
    # _axis1_plt4 = _retro_df.plot(y=_colname, kind='kde', secondary_y=True, mark_right=False, legend=False, lw=0.4,color='magenta', ax=_axis_plt4)

    _axis_plt4.hist(x=_retro_df[_colname], bins=_bins, histtype='bar', color='lightgray')

    _axis1_plt4 = _axis_plt4.twinx()
    #_retro_df.plot(y=_colname, kind='kde', legend=False, ax=_axis1_plt4, color='dimgray', mark_right=True, lw=0.7, alpha=0.8)
    sns.distplot(_retro_df[_colname], hist=False, bins=_bins, ax=_axis1_plt4, color='dimgray', kde_kws={'lw': 0.7}, hist_kws={'alpha': 0.8})


    _axis_plt4.set_xlim(0, 105)

    _axis_plt4.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(_bins[0::5], nbins=21))
    _axis_plt4.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))

    _axis_plt4.tick_params(axis='x', which='both', length=1, width=0.5, labelbottom=True, bottom=True, labelsize=3, direction='out', pad=2)
    _axis_plt4.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=4, labelleft=True, left=True, direction='out', pad=2)

    _axis_plt4.set_title(_plot_label, fontsize=4.4)


    _axis_plt4.set_xlabel('% Mapped / Aligned Reads', labelpad=1, fontsize=3)

    if _current_sample <= _exonMapping_cutoff_fail:
        print("Sample flagged for FAILURE!!!")
        insert_flag_fail(_axis_plt4)

    elif (_current_sample > _exonMapping_cutoff_fail) and (_current_sample <= _exonMapping_cutoff_warn):
        print("Sample flagged for WARNING!!!")
        insert_flag_warn(_axis_plt4)

    else:
        print("Sample PASSED cutoffs!")
        pass

    ### Adding cutoff markers
    _axis_plt4.plot(_exonMapping_cutoff_fail, _axis_plt4.get_ylim()[1] - (_axis_plt4.get_ylim()[1] / 10), marker='v', ms=0.8, c='red')
    _axis_plt4.text(_exonMapping_cutoff_fail, _axis_plt4.get_ylim()[1] - (_axis_plt4.get_ylim()[1] / 20), 'Fail', fontsize=2, color='red', horizontalalignment='center')

    _axis_plt4.plot(_exonMapping_cutoff_warn, _axis_plt4.get_ylim()[1] - (_axis_plt4.get_ylim()[1] / 10), marker='v', ms=0.8, c='gold')
    _axis_plt4.text(_exonMapping_cutoff_warn, _axis_plt4.get_ylim()[1] - (_axis_plt4.get_ylim()[1] / 20), 'Warn', fontsize=2, color='gold', horizontalalignment='center')

    _axis_plt4.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    _axis_plt4.set_ylabel('Frequency', labelpad=2, fontsize=3)

    # SECONDARY AXIS FOR KDE
    _axis1_plt4.get_yaxis().set_ticks([])
    _axis1_plt4.yaxis.label.set_visible(False)
    # _axis1_plt4.set_ylabel('Kernel Density Estimate', fontsize=4, labelpad=3)

    _axis_plt4.spines['top'].set_visible(False)
    _axis_plt4.spines['right'].set_visible(False)
    _axis_plt4.spines['left'].set_visible(True)
    _axis_plt4.spines['bottom'].set_visible(True)
    _axis_plt4.spines['left'].set_color('black')
    _axis_plt4.spines['bottom'].set_color('black')
    _axis_plt4.spines['bottom'].set_linewidth(0.55)
    _axis_plt4.spines['left'].set_linewidth(0.55)

    _axis1_plt4.spines['top'].set_visible(False)
    _axis1_plt4.spines['right'].set_visible(False)
    _axis1_plt4.spines['left'].set_visible(False)
    _axis1_plt4.spines['bottom'].set_visible(False)

    # _axis_plt4.grid(b=False)
    # _axis1_plt4.grid(b=False)

    if _current_sample > _lib_mean:
        _adj_flag = True
        _rotation_current = 270
        _rotation_mean = 90
    elif _current_sample < _lib_mean:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270
    else:
        _adj_flag = False
        _rotation_current = 90
        _rotation_mean = 270

    # Current Sample Line and Label
    _line1 = _axis_plt4.axvline(x=_current_sample, alpha=0.8, color='red', linestyle='-', linewidth=0.5, label='{:2.2f}%'.format(_current_sample))

    if _adj_flag:
        _axis_plt4.text(_current_sample + 0.5, (_axis_plt4.get_ylim()[1] / 2), '{:2.2f}%'.format(_current_sample), rotation=_rotation_current, fontsize=2.5, zorder=2)
    else:
        _axis_plt4.text(_current_sample - 2, (_axis_plt4.get_ylim()[1] / 2), '{:2.2f}%'.format(_current_sample), rotation=_rotation_current, fontsize=2.5, zorder=2)

    # Current Library Mean Line and Label
    _line2 = _axis_plt4.axvline(x=_lib_mean, alpha=0.8, color='indigo', linestyle='--', linewidth=0.5, label='{:2.2f}%'.format(_lib_mean))

    if _adj_flag:
        _axis_plt4.text(_lib_mean - 2, ((_axis_plt4.get_ylim()[1] / 2) + 1), '{:2.2f}%'.format(_lib_mean), rotation=_rotation_mean, fontsize=2.5, zorder=2)
    else:
        _axis_plt4.text(_lib_mean + 0.5, ((_axis_plt4.get_ylim()[1] / 2) + 1), '{:2.2f}%'.format(_lib_mean), rotation=_rotation_mean, fontsize=2.5, zorder=2)

    ## Superimpose the Kernel Density Estimate line over the distribution
    _kde_line = matplotlib.lines.Line2D([0], [0], color="dimgray", linewidth=0.5, linestyle='-')

    _axis_plt4.legend([_line1, _line2], ["Current Sample", "Batch Mean"], loc='best', frameon=False, fontsize=3, ncol=1)
    _axis_plt4.set_facecolor('white')

    return _figure



#### Plot 5: rRNA Scatter ####
def plotScatter_rRNA(_in_tup, _userDf, _background_df, _pos, _f=None):
    print(_in_tup)

    # print(min(_background_df['Num_Uniquely_Aligned_rRNA']))
    # print(_background_df['Num_Uniquely_Aligned'])

    # print(_in_tup)
    # print(_in_tup[1])#Sample Name from Input Tuple
    # print(_userDf)

    _plotter_df = pd.concat([_background_df, _userDf])

    # Assign color for current project's library (all samples in the current project)
    _plotter_df["scatter_color"] = np.where(_plotter_df["Sample"].isin(_userDf["Sample"]), "indigo", "darkgray")

    print(_plotter_df)

    # print(_background_df.loc[_background_df["Sample"] == _in_tup[1], "Sample"])

    # Assign separate color for current sample on each page
    _plotter_df.loc[_plotter_df["Sample"] == _in_tup[1], "scatter_color"] = 'r'
    print(_plotter_df)

    # print(_background_df['Num_Uniquely_Aligned'])
    # print(_background_df['Num_Uniquely_Aligned_rRNA'])

    if not _f is None:
        plt.gcf()

    _ax = plt.subplot(3, 2, _pos)

    _ax.scatter(x=_plotter_df['Num_Uniquely_Aligned'], y=_plotter_df['Num_Uniquely_Aligned_rRNA'], s=0.5, c=_plotter_df["scatter_color"])

    ## Regression line (gradient slope)

    X = _plotter_df.loc[:, "Num_Uniquely_Aligned"].values.reshape(-1, 1)
    Y = _plotter_df.loc[:, "Num_Uniquely_Aligned_rRNA"].values.reshape(-1, 1)
    linear_regressor = LinearRegression()
    linear_regressor.fit(X, Y)
    Y_pred = linear_regressor.predict(X)

    _ax.plot(X, Y_pred, c='dimgray', linewidth=0.7, linestyle='-', alpha=1)

    _ax.set_title("Ribosomal RNA", fontsize=4.4)

    _ax.tick_params(axis='x', which='both', length=1, width=0.5, labelrotation=0, labelbottom=True, bottom=True, labelsize=3, direction='out', pad=1)
    _ax.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=3, labelleft=True, left=True, direction='out', pad=2)

    _ax.set_xlabel("Uniquely Aligned Reads (Millions)", fontsize=4, labelpad=2)
    _ax.set_ylabel("rRNA Reads (Millions)", fontsize=4, labelpad=2)

    # _ax.set_xlim(0, max(_background_df['Num_Uniquely_Aligned_rRNA']) + 1000000)
    # _ax.set_ylim(0, max(_background_df['Num_Uniquely_Aligned']) + 1000000)
    # Set limits on X-axis
    # _ax.set_xlim([0, max(_background_df['Num_Uniquely_Aligned_rRNA'])+2000000])
    # _ax.axis('scaled')

    x_bottom, x_top = plt.xlim()
    print(x_bottom, x_top)

    y_bottom, y_top = plt.ylim()
    print(y_bottom, y_top)

    print(_in_tup)

    # _ax.axis('image')

    # Set limits on X-axis
    _ax.set_xlim([x_bottom - 400000, x_top + 3000000], emit=True)

    # Set limits on Y-axis
    _ax.set_ylim([y_bottom, y_top + 7500000], emit=True)

    # Plotting the ratio line (FAIL_slope=0.50 (50%), WARN_slope=0.35 (35%))
    _slope_warn = 0.35
    _slope_fail = 0.50

    _slope_current = float(_in_tup[7] / _in_tup[4])
    print(_slope_current)

    if _slope_current >= _riboScatter_cutoff_fail:
        print("Sample FAILED!")
        insert_flag_fail(_ax)
    elif _slope_current >= _riboScatter_cutoff_warn:
        print("WARNING!!!")
        insert_flag_warn(_ax)
    else:
        print("Sample PASSED.")

    xmin, xmax = _ax.get_xlim()
    line_x0 = 0
    line_x1 = xmax

    line_y0 = 0
    line_y1_warn = _slope_warn * (line_x1 - line_x0) + line_y0
    line_y1_fail = _slope_fail * (line_x1 - line_x0) + line_y0

    _ax.plot([line_x0, line_x1], [line_y0, line_y1_warn], c='orange', linewidth=0.3, linestyle='--', alpha=0.3, label="Warn")
    _ax.plot([line_x0, line_x1], [line_y0, line_y1_fail], c='r', linewidth=0.3, linestyle='--', alpha=0.3, label="Fail")

    # _ax.annotate('Warn', xy=(45000000, line_y1_warn), xytext=(45000000, 15000000), fontsize=3, color='orange', ha='center', va='center')
    _ax.annotate('Warn', xy=(line_x1 - 1000000, line_y1_warn - 1000000), fontsize=3, color='orange', ha='right', va='center', rotation=13.5)
    # arrowprops=dict(arrowstyle='<->', connectionstyle='arc3,rad=0', lw=0.3, ls='-'),

    # _ax.annotate('Fail', xy=(45000000, line_y1_fail), xytext=(45000000, 22500000), fontsize=3, color='r', ha='center', va='center')
    _ax.annotate('Fail', xy=(line_x1 - 3000000, line_y1_fail - 2000000), fontsize=3, color='r', ha='right', va='center', rotation=17)
    # arrowprops=dict(arrowstyle='<-', connectionstyle='arc3,rad=0', lw=0.4, ls='-'),

    ## PASS deprecated
    # _ax.annotate('Pass', xy=(45000000, line_y1_warn), xytext=(45000000, 7500000), fontsize=3, color='green', ha='center', va='center')
    # arrowprops=dict(arrowstyle='-', connectionstyle='arc3,rad=0', lw=0.4, ls='-'),

    # Set axes margins for padding on both axes
    _ax.margins(0.01)

    # _ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=False))
    _ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2500000))
    # _ax.xaxis.set_major_locator(matplotlib.ticker.IndexLocator(base=2500000, offset=0))
    _ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_scatter_million))

    _ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2500000))
    # _ax.yaxis.set_major_locator(matplotlib.ticker.IndexLocator(base=2500000, offset=0))
    _ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_scatter_million))

    _ax.spines['top'].set_visible(False)
    _ax.spines['bottom'].set_visible(True)
    _ax.spines['bottom'].set_color('black')
    _ax.spines['right'].set_visible(False)

    _ax.spines['left'].set_visible(True)
    _ax.spines['left'].set_color('black')

    _ax.spines['left'].set_linewidth(0.55)
    _ax.spines['bottom'].set_linewidth(0.55)
    _ax.set_facecolor('white')

    print(_ax.get_data_ratio())
    # _ax.set_aspect(1/_ax.get_data_ratio()*0.40)
    _ax.set_aspect('auto', adjustable='box', anchor='SW')

    _curr_samp = matplotlib.lines.Line2D([0], [0], color='w', markerfacecolor='r', marker='o', linewidth=1, markersize=3.5)
    _curr_lib = matplotlib.lines.Line2D([0], [0], color='w', markerfacecolor='indigo', marker='o', linewidth=1, markersize=3.5)
    _historic_data = matplotlib.lines.Line2D([0], [0], color='w', markerfacecolor='darkgray', marker='o', linewidth=1, markersize=3.5)
    _regression_gradient = matplotlib.lines.Line2D([0], [0], color='black', linewidth=0.6)

    _ax.legend([_curr_samp, _curr_lib], ['Current Sample', 'Batch Samples'], loc='best', frameon=False, ncol=1, fontsize=3)

    return _f


#### Plot 6: Sequence Contamination - Violin Plot ####
def plotViolin_dualAxis(_input_tup, _userDf, _background_df, _position, _f=None):
    
    
    ## Load the sequence contamination levels for the background data
    _contaminant_df_untrim = _background_df[["Percent_Overrepresented_Seq_Untrimmed", "Percent_Adapter_Content_Untrimmed"]]
    _contaminant_df_trim = _background_df[["Percent_Overrepresented_Seq_Trimmed", "Percent_Adapter_Content_Trimmed"]]


    _contaminant_df_untrim.columns = ["Overrepresented", "Adapter"]
    _contaminant_df_trim.columns = ["Overrepresented", "Adapter"]


    _contaminant_melt_untrim = pd.melt(_contaminant_df_untrim, var_name="Contamination_Metric", value_name="Percent")
    _contaminant_melt_trim = pd.melt(_contaminant_df_trim, var_name="Contamination_Metric", value_name="Percent")


    print(_input_tup)


    _current_overrep_untrim = _input_tup[8]
    _current_adapter_untrim = _input_tup[9]

    _current_overrep_trim = _input_tup[10]
    _current_adapter_trim = _input_tup[11]


    ## Remove the current batch mean from the USER dataframe
    _user_minusBatchMean_df = _userDf.drop(_userDf.tail(1).index)
    

    _mean_overrep_untrim = _user_minusBatchMean_df.loc[:, 'Percent_Overrepresented_Seq_Untrimmed'].mean()
    _mean_adapter_untrim = _user_minusBatchMean_df.loc[:, 'Percent_Adapter_Content_Untrimmed'].mean()

    _mean_overrep_trim = _user_minusBatchMean_df.loc[:, 'Percent_Overrepresented_Seq_Trimmed'].mean()
    _mean_adapter_trim = _user_minusBatchMean_df.loc[:, 'Percent_Adapter_Content_Trimmed'].mean()



    ## Define color palette
    _contaminant_pal = {"Overrepresented": "lightgray", "Adapter": "gray"}

    if not _f is None:
        plt.gcf()


    _gridsp = matplotlib.gridspec.GridSpec(6, 2, figure=_f)


    _axis = _f.add_subplot(_gridsp[4, 1:])
    _axis2 = _f.add_subplot(_gridsp[5, 1:])



    sns.violinplot(x="Percent", y="Contamination_Metric", data=_contaminant_melt_untrim, palette=_contaminant_pal,
                   inner=None,
                   ax=_axis, linewidth=0.3, orient="h", scale="count")

    # _axis2 = _axis.twiny()

    sns.violinplot(x="Percent", y="Contamination_Metric", data=_contaminant_melt_trim, palette=_contaminant_pal,
                   inner=None,
                   ax=_axis2, linewidth=0.3, orient="h", scale="count")

    _axis.set_title("Sequence Contamination", fontsize=4.4, pad=0)

    _axis.tick_params(axis='x', which='both', length=1, width=0.5, labelsize=3, labelbottom=True, bottom=True, direction='out', pad=2)
    _axis.tick_params(axis='y', which='both', length=1, width=0.5, labelrotation=0, labelleft=True, left=True, labelsize=2.5, direction='out', pad=1)

    _axis2.tick_params(axis='x', which='both', length=1, width=0.5, labelsize=3, labelbottom=True, bottom=True, direction='out', pad=2)
    _axis2.tick_params(axis='y', which='both', length=1, width=0.5, labelrotation=0, labelleft=True, left=True, labelsize=2.5, direction='out', pad=1)

    _axis.set_xlabel("Pre-trim (%)", fontsize=3.5, labelpad=0.5)
    _axis.set_ylabel("")

    _axis2.set_xlabel("Post-trim (%)", fontsize=3.5, labelpad=0.5)
    _axis2.set_ylabel("")

    # _axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_contaminant))
    print(_axis.get_xlim()[0], _axis.get_xlim()[1])
    _axis.set_xlim(_axis.get_xlim()[0] + 1, _axis.get_xlim()[1] + 5)

    _axis2.set_xlim(_axis2.get_xlim()[0], _axis2.get_xlim()[1] + 2)

    # print(_axis.get_xticks(), _axis.get_xticklabels())

    _axis.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    _axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))

    _axis2.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    _axis2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))

    _axis.set_yticklabels(['Overrepresented', 'Adapter'])
    _axis2.set_yticklabels(['Overrepresented', 'Adapter'])

    _axis.spines['top'].set_visible(False)
    _axis.spines['bottom'].set_visible(True)
    _axis.spines['bottom'].set_color('black')
    _axis.spines['right'].set_visible(False)

    _axis.spines['left'].set_visible(True)
    _axis.spines['left'].set_color('black')

    _axis.spines['left'].set_linewidth(0.55)
    _axis.spines['bottom'].set_linewidth(0.55)

    _axis2.spines['top'].set_visible(False)
    _axis2.spines['bottom'].set_visible(True)
    _axis2.spines['bottom'].set_color('black')
    _axis2.spines['right'].set_visible(False)

    _axis2.spines['left'].set_visible(True)
    _axis2.spines['left'].set_color('black')

    _axis2.spines['left'].set_linewidth(0.55)
    _axis2.spines['bottom'].set_linewidth(0.55)

    _x_bottom, _x_top = _axis.get_xlim()
    print(_x_bottom, _x_top)
    print(_axis.get_xlim())

    _y_bottom, _y_top = _axis.get_ylim()
    print(_y_bottom, _y_top)

    if (_current_overrep_trim or _current_adapter_trim) >= _violin_cutoff_fail:
        print("FAILED!")
        insert_flag_fail(_axis)
    elif (_current_overrep_trim or _current_adapter_trim) >= _violin_cutoff_warn:
        print("WARNING!!!")
        insert_flag_warn(_axis)
    else:
        print("PASSED!!")
        pass

    '''
    
    ### For separating the cutoffs for Overrepresented and Adapter ###
    if _current_overrep_trim >= _violin_cutoff_overrep_fail:
        print("Overrepresented FAILED!")
        insert_flag_fail(_axis)
    elif _current_overrep_trim >= _violin_cutoff_overrep_warn:
        print("Overrepresented WARNING!!!")
        insert_flag_warn(_axis)
    else:
        print("PASSED!!")
        pass

    if _current_adapter_trim >= _violin_cutoff_adapter_fail:
        print("Adapter FAILED!")
        insert_flag_fail(_axis)
    elif _current_adapter_trim >= _violin_cutoff_adapter_warn:
        print("Adapter WARNING!!!")
        insert_flag_warn(_axis)
    else:
        print("PASSED!!")
        pass
    '''


    ### Adding cutoff markers
    _axis3 = _axis2.twinx()

    # _axis3.xaxis.set_ticks([])
    _axis3.yaxis.set_ticks([])

    _axis3.xaxis.label.set_visible(False)
    _axis3.yaxis.label.set_visible(False)

    _axis3.set_ylim(_axis2.get_ylim()[0], _axis2.get_ylim()[1])

    _axis3.plot(_violin_cutoff_fail, -0.7, marker='v', ms=1, c='red', clip_on=False)
    _axis3.text(_violin_cutoff_fail, -0.88, 'Fail', fontsize=3, color='red', horizontalalignment='center')

    _axis3.plot(_violin_cutoff_warn, -0.7, marker='v', ms=1, c='gold', clip_on=False)
    _axis3.text(_violin_cutoff_warn, -0.88, 'Warn', fontsize=3, color='gold', horizontalalignment='center')

    _axis3.spines['top'].set_visible(False)
    _axis3.spines['right'].set_visible(False)
    _axis3.spines['bottom'].set_visible(False)
    _axis3.spines['left'].set_visible(False)



    _line_overrep_untrim = _axis.axvline(x=_current_overrep_untrim, ymin=0.55, ymax=0.95, alpha=0.8, color='red',
                                         linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_overrep_untrim))
    _line_mean_overrep_untrim = _axis.axvline(x=_mean_overrep_untrim, ymin=0.55, ymax=0.95, alpha=0.8, color='indigo',
                                              linestyle='--', linewidth=0.35,
                                              label='{:.2f}%'.format(_mean_overrep_untrim))

    _line_adapter_untrim = _axis.axvline(x=_current_adapter_untrim, ymin=0.05, ymax=0.45, alpha=0.8, color='red',
                                         linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_adapter_untrim))
    _line_mean_adapter_untrim = _axis.axvline(x=_mean_adapter_untrim, ymin=0.05, ymax=0.45, alpha=0.8, color='indigo',
                                              linestyle='--', linewidth=0.35,
                                              label='{:.2f}%'.format(_mean_adapter_untrim))

    _line_overrep_trim = _axis2.axvline(x=_current_overrep_trim, ymin=0.55, ymax=0.95, alpha=0.8, color='red',
                                        linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_overrep_trim))
    _line_mean_overrep_trim = _axis2.axvline(x=_mean_overrep_trim, ymin=0.55, ymax=0.95, alpha=0.8, color='indigo',
                                             linestyle='--', linewidth=0.35, label='{:.2f}%'.format(_mean_overrep_trim))

    _line_adapter_trim = _axis2.axvline(x=_current_adapter_trim, ymin=0.05, ymax=0.45, alpha=0.8, color='red',
                                        linestyle='-', linewidth=0.35, label='{:.2f}%'.format(_current_adapter_trim))
    _line_mean_adapter_trim = _axis2.axvline(x=_mean_adapter_trim, ymin=0.05, ymax=0.45, alpha=0.8, color='indigo',
                                             linestyle='--', linewidth=0.35, label='{:.2f}%'.format(_mean_adapter_trim))

    ## Add dividing line between pre-trim and post-trim along with descriptive text
    # _dividing_line = _axis.axhline(y=1.5, xmin=0, xmax=1, color="darkslategray", alpha=0.3, linestyle='--', linewidth=0.3)

    # _axis.annotate('Pre-Trim', xy=(55, 1.3), fontsize=2.5, color='r', ha='center', va='center', rotation=0)
    # _axis.annotate('Post-Trim', xy=(55, 3.3), fontsize=2.5, color='r', ha='center', va='center', rotation=0)

    # _axis.text(0.9, 0.5, 'Pre-Trim', fontsize=2.5, zorder=2)
    # _axis.text(0.9, 0.2, 'Post-Trim', fontsize=2.5, zorder=2)

    _axis.legend([_line_overrep_trim, _line_mean_overrep_trim], ["Current Sample", "Batch Mean"], loc='best',
                 frameon=False, ncol=1, fontsize=3)

    _axis.set_facecolor('white')
    _axis2.set_facecolor('white')
    
    plt.subplots_adjust(hspace=0)
    

    return _f






#### Plot 7 : GeneBody Coverage Plot

def plotGC(_ipTuple, _coverage_df, _position, _plot_title, _fig=None):
    print(_ipTuple[1])

    if not _fig is None:
        plt.gcf()

    _axis = _fig.add_subplot(4, 2, _position)


    # Calculate mean GeneBody Coverage for the entire library
    _mean_df = pd.DataFrame()
    _mean_df['gc_mean'] = _coverage_df.mean(axis=1)
    # print(_mean_df)

    # Statistical Tests: Wilcoxon signed rank test, Kolmogorov-Smirnov 2-sample test and Pearson Correlation
    _wilcox_stat, _wilcox_pval = stats.wilcoxon(_coverage_df[_ipTuple[1]], _mean_df['gc_mean'], correction=True)
    print(_wilcox_stat, _wilcox_pval)

    ## KS-2sample test
    _ks_stat, _ks_pval = stats.ks_2samp(_coverage_df[_ipTuple[1]], _mean_df['gc_mean'])
    print("KS-2samp test")
    print(_ks_stat, _ks_pval)

    _pearson_corr = _coverage_df[_ipTuple[1]].corr(_mean_df['gc_mean'])
    print("PEARSON CORR :", _pearson_corr)

    ### Calculate Confidence Interval for the mean GC line
    # _sigma = stats.t.interval(0.05, len(_mean_df['gc_mean'])-1, loc=_mean_df['gc_mean'], scale=stats.sem(_mean_df['gc_mean']))
    _err = stats.sem(_mean_df['gc_mean']) * stats.t.ppf((1 + 0.95) / 2, len(_mean_df) - 1)
    print(_err)

    # Plot current sample with library mean
    _x = np.arange(1, 101, 1)

    _axis.plot(_x, _coverage_df[_ipTuple[1]], color='red', linewidth=0.5, linestyle='-', alpha=0.8)
    _axis.plot(_x, _mean_df['gc_mean'], color='indigo', linewidth=0.5, linestyle='--', alpha=0.8)

    _axis.fill_between(_x, _mean_df['gc_mean'] - _err, _mean_df['gc_mean'] + _err, facecolor='yellow', alpha=0.5)

    _axis.tick_params(axis='x', which='both', length=1, width=0.5, labelbottom=True, bottom=True, labelsize=3,
                      direction='out', pad=2)
    _axis.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=3, labelleft=True, left=True,
                      direction='out', pad=2)

    _axis.set_xlim(0, 105)
    _axis.set_title(_plot_title, fontsize=4.4)

    _uni_arrw = u"\u2192"
    _axis.set_xlabel("Gene Percentile (5' " + _uni_arrw + " 3')", fontsize=3, labelpad=2)
    _axis.set_ylabel("Coverage", fontsize=3, labelpad=2)

    _axis.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(np.arange(0, 101, 1)[0::5]))
    _axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))
    # _axis.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    _axis.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt_cov))

    _axis.spines['top'].set_visible(False)
    _axis.spines['right'].set_visible(False)
    _axis.spines['left'].set_visible(True)
    _axis.spines['bottom'].set_visible(True)
    _axis.spines['left'].set_color('black')
    _axis.spines['bottom'].set_color('black')
    _axis.spines['bottom'].set_linewidth(0.55)
    _axis.spines['left'].set_linewidth(0.55)

    ## KS-2sample flagging
    if _ks_pval <= 0.05:
        print("KS-2samp Test : FAILED!\n")
        insert_flag_fail(_axis)
    elif _ks_pval <= 0.1:
        print("KS-2samp Test : WARNING!\n")
        insert_flag_warn(_axis)
    else:
        print("KS-2samp Test : PASSED!\n")
        pass

    ### Wilcoxon flagging
    # if _wilcox_pval <= 0.05:
    #    print("Wilcoxon Test : FAILED!\n")
    #    insert_flag_fail(_axis)
    # elif _wilcox_pval <= 0.1:
    #    print("Wilcoxon Test: WARNING!\n")
    #    insert_flag_warn(_axis)
    # else:
    #    print("Wilcoxon Test : PASSED!\n")
    #    pass

    _current_sample_line = matplotlib.lines.Line2D([0], [0], color="red", linewidth=0.5, linestyle='-', alpha=0.8)
    _library_line = matplotlib.lines.Line2D([0], [0], color="indigo", linewidth=0.5, linestyle='--', alpha=0.8)

    _extra_confidenceInterval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='yellow', fill=True,
                                                             edgecolor='yellow', linewidth=1.2, alpha=0.5)
    _extra_ksPval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None', linewidth=0)
    # _extra_wilcox = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None', linewidth=0)

    _axis.legend([_current_sample_line, _library_line, _extra_confidenceInterval, _extra_ksPval],
                 ["Current Sample", "Library Mean",
                  "95% Confidence Interval", "KS-2sample Pvalue: " + str(round(_ks_pval, 5))],
                 loc='best', frameon=False, fontsize=3, ncol=1)

    return _fig



#### Plot 8 : Gene Expression Distribution Plot 
def plotNegBin(_ipTuple, _hist_df, _user_df, _pos, _plot_title, _f=None):
    _index_array = _hist_df.iloc[:, 0]
    print(_index_array)

    _low_vals = []
    _high_vals = []

    for _i in _index_array:
        # print(_i.strip('(').strip(']').split(','))
        _low_vals.append(float(_i.strip('(').strip(']').split(',')[0]))
        _high_vals.append(float(_i.strip('(').strip(']').split(',')[1]))

    print(_low_vals)

    _col_index = pd.IntervalIndex.from_arrays(_low_vals, _high_vals, closed='right')
    print(len(_col_index))

    _col_index = _col_index[:25]

    print(_col_index)

    # _x_vals = np.arange(0, len(_col_index) / 2, 0.5)
    _x_vals = np.arange(0, len(_col_index) / 2, 0.5)
    print(_x_vals)
    print(len(_x_vals))

    ## Preparing the data_df and libMean_df for all bins
    _data_df = _hist_df.drop(['Unnamed: 0'], axis=1)
    _libMean_df = pd.DataFrame()
    _libMean_df['Mean'] = _data_df.iloc[:, :-1].mean(numeric_only=True, axis=1)
    print(_libMean_df)

    ## Dropping the last 8 bins from data_df and libMean_df to get focused resolution
    _data_df_dropped = _data_df.drop(_data_df.tail(8).index)
    _mean_df_dropped = _libMean_df.drop(_libMean_df.tail(8).index)

    _max_df = pd.DataFrame()

    # print(_mean_df)
    # _mean_df['hist_mean'] = _data_df.mean(axis=0)

    # print(_mean_df)
    # _mean_df.drop(_mean_df.tail(8).index, inplace=True)

    _max_df['max_val'] = _data_df.max(axis=1)
    print(_max_df['max_val'].max())

    ### Statistical Tests
    # Wilcoxon signed rank test
    # _wilcox_stat, _wilcox_pval = stats.wilcoxon(_data_df[_ipTuple[1]], _mean_df['hist_mean'], correction=True)
    # print(_wilcox_stat, _wilcox_pval)

    # Kolmogorov-Smirnov 2-sample test for comparing the current distribution to the mean of the library
    # _ksTest_stat, _ksTest_pval = stats.ks_2samp(_data_df[_ipTuple[1]], _mean_df['hist_mean'])
    # print("KS 2sample Test :")
    # print(_ksTest_stat, _ksTest_pval)

    ## Expressed genes correlation test Pearson's Correlation
    _mean_array = _mean_df_dropped.Mean.values
    _current_samp_array = _data_df_dropped[_ipTuple[1]].values

    print(_mean_array)
    print(_current_samp_array)

    ## Zscore calculated for each sample in the dataframe
    # _zscore_data_df = _data_df_dropped.apply(stats.zscore)
    # print(_zscore_data_df)

    ## ZTest for mean and current sample against the normal distribution to get pvalue
    _ztest_stat_raw, _ztest_pval_raw = ztest_prob(_current_samp_array, _mean_array, 0)

    # _ztest_stat_1samp, _ztest_pval_1samp = ztest_prob(_current_samp_dist, None, 0)

    # _ztest_pval_zscore = ztest_prob(_zscore_data_df[_ipTuple[1]].values, _zscore_mean_df.hist_mean.values)
    # print(_zscore_mean_df)

    # print(_zscore_data_df.loc[:,_ipTuple[1]])
    # print(_zscore_mean_df)
    # print(_current_samp_dist)

    # _pearson_corr = _data_df[_ipTuple[1]].corr(_mean_df['hist_mean'])
    # _pearson_corr, _pearson_pval = stats.stats.pearsonr(_current_samp_hist, _mean_hist)
    # print("Pearson P-value : ", _pearson_pval)

    if not _f is None:
        plt.gcf()

    _ax = _f.add_subplot(4, 2, _pos)

    _col_names = [_cl for _cl in _data_df_dropped.columns]
    print(_col_names)

    ## Plotting all current library distributions with the current sample highlighted on each page
    for _col in _col_names:

        if _col == _ipTuple[1]:
            plt.plot(_x_vals, _data_df_dropped[_col], color='red', linewidth=0.5, linestyle='-', zorder=24)
        else:
            plt.plot(_x_vals, _data_df_dropped[_col], color='silver', linewidth=0.5, linestyle='-')

    # Removed plotting current line separately since the above loop
    # _ax.plot(_x_vals, _data_df_dropped[_ipTuple[1]], color='red', linewidth=0.5, linestyle='-', alpha=0.8)

    ## Plotting the mean distribution
    _ax.plot(_x_vals, _mean_df_dropped['Mean'], color='indigo', linewidth=0.5, linestyle='--', alpha=0.8, zorder=23)

    _ax.tick_params(axis='x', which='both', length=1, width=0.5, labelrotation=60, labelbottom=True, bottom=True,
                    labelsize=2, direction='out', pad=1)
    _ax.tick_params(axis='y', which='both', length=1, width=0.5, labelsize=3, labelleft=True, left=True,
                    direction='out', pad=2)

    # _ax.set_xlim(0, (len(_x_vals) / 2) + 1)
    _ax.set_xlim(0, 14)

    _ax.set_ylim(0, 4000)
    # _ax.set_ylim(0, _max_df['max_val'].max() + 300)

    _ax.set_title(_plot_title, fontsize=4.4)

    _ax.set_xlabel("Gene Expression Value (log2(CPM)+1)", fontsize=3, labelpad=2)
    _ax.set_ylabel("Frequency", fontsize=3, labelpad=2)

    _ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(_x_vals))
    _ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(_col_index))
    _ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(500))
    _ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(fmt))

    # _ax.xticklabels(str(_col_index), fontsize=3)

    _ax.spines['top'].set_visible(False)
    _ax.spines['right'].set_visible(False)
    _ax.spines['left'].set_visible(True)
    _ax.spines['bottom'].set_visible(True)

    _ax.spines['left'].set_color('black')
    _ax.spines['left'].set_linewidth(0.55)
    # _ax.spines['left'].set_smart_bounds(True)

    _ax.spines['bottom'].set_color('black')
    _ax.spines['bottom'].set_linewidth(0.55)
    # _ax.spines['bottom'].set_smart_bounds(True)

    ## Flagging based on Ztest
    if _ztest_pval_raw <= 0.05:
        print("ZTest : FAILED!")
        insert_flag_fail(_ax)
    elif _ztest_pval_raw <= 0.1:
        print("ZTest : WARNING!")
        insert_flag_warn(_ax)
    else:
        pass

    _current_samp_line = matplotlib.lines.Line2D([0], [0], color="red", linewidth=0.5, linestyle='-', alpha=0.8)
    _lib_line = matplotlib.lines.Line2D([0], [0], color="indigo", linewidth=0.5, linestyle='--', alpha=0.8)

    # _extra_wilcox_stat = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None',linewidth=0)
    # _extra_wilcox_pval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None', linewidth=0)
    # _extra_ksStat = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None', linewidth=0)
    # _extra_ksPval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None', linewidth=0)

    _extra_Ztest_stat = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None',
                                                     linewidth=0)
    _extra_Ztest_Pval = matplotlib.patches.Rectangle((0, 0), 1, 1, facecolor='w', fill=False, edgecolor='None',
                                                     linewidth=0)

    _ax.legend([_current_samp_line, _lib_line, _extra_Ztest_Pval],
               ["Current Sample", "Library Mean", "ZTest Pvalue : " + str(round(_ztest_pval_raw, 5))], loc='best',
               frameon=False, fontsize=3, ncol=1)

    return _f






if __name__ == "__main__":
    main()