import matplotlib
import matplotlib.gridspec as gridspec
#matplotlib.use('Agg')
import seaborn
import matplotlib.pyplot as plt
#import sqlite3
import os
import sys
import json
import seaborn as sns
import time
import shutil
import glob
#import sqlalchemy as sq
import numpy as np
import pandas as pd
import xlsxwriter
from itertools import zip_longest
from scipy.stats import norm
from sklearn.preprocessing import LabelEncoder
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.preprocessing import MinMaxScaler



####### Add all the data loading, cleaning and other helper functions here ##########

####### Add all the plotting functions here for all panels ########





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