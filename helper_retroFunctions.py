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



#### Plot 1: Input Size #### (Restructure the names of the dataframes)
def plotHist_ipSize(_in_tuple, _userDf, _background_df, _pos, _f=None):
    # print(_background_df.Input_Size)
    bin = np.arange(0, 100 + 1, 1)

    print(_background_df)
    print(_background_df.columns)

    ### Add a hard-clip limit of 60 million reads to the master['Input_Size']
    _background_df.loc[:, 'Input_Size'] = _background_df.loc[:, 'Input_Size'].clip(upper=40)
    print(_background_df.Input_Size)

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

    _ax = _f.add_subplot(4, 2, _pos)

    # _ax = _background_df['Input_Size'].plot(kind='hist', bins=_bins, ax=plt.gca())
    # _ax1 = _background_df.plot(y='Input_Size', kind='kde', secondary_y=True, mark_right=True, legend=False, lw=0.5, ax=_ax, color='magenta')

    # _ax.hist(_background_df['Input_Size'], bins=_bins)
    _background_df["Input_Size"].plot(kind='hist', bins=_bins, ax=_ax, color='lightgray')

    _ax1 = _ax.twinx()
    _background_df.plot(y="Input_Size", kind='kde', legend=False, ax=_ax1, color='dimgray', mark_right=True, lw=0.7,
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



if __name__ == "__main__":
    main()