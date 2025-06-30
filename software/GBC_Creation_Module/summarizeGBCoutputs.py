'''
Author:     Tyler Therron
Maintainer: Samuel Hamilton
Email:      samuelhamilton2024@u.northwestern.edu
'''

import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') # set backend to agg
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.preprocessing import normalize
import argparse
from sklearn.preprocessing import MinMaxScaler
from matplotlib.backends.backend_pdf import PdfPages
import seaborn

def geneBody_summaryTable(_input_txt_dir, _out_dir, _project_name):
# define the directory path
    my_path = _input_txt_dir
 
    # get the list of files in the directory
    files_list = os.listdir(my_path) 

    # create output directory for figure
   #  - no need to make new directory
    # definer scaler as object
    scaler = MinMaxScaler()   
# create an empty list 
    df_list = []

    df_list2 = []
#empty pandas dataframe to append to
    genebody_df = pd.DataFrame()

# iterate through each file in the directory
    for file in files_list:
    # read the data of each file
        with open(my_path + "/" + file, 'r') as f: # 2/6/23 - added slash because of file location error
            data = pd.read_csv(f, sep="\t")
        # read in the tab delimited data into the 'data' object and make a dataframe from it
            pd.DataFrame(data)
        
        # append data to list
            df_list.append(data)
        
    result = pd.concat(df_list, axis=0, join='outer', ignore_index=True)

# Rearranging index
    result.index = np.arange(1, len(result) + 1)

# bam files have the '_chr' cut from the name
    result['Percentile'] = result['Percentile'].str.replace('_chr', '')
    
    # this is where the changes begin with the new normalization methods after talking with Debbie
    result_transform = result
    result_transform = result_transform.drop('Percentile',axis=1)

    result_scaled=pd.DataFrame(scaler.fit_transform(result_transform.T).T,columns=result_transform.columns)
    
    # MaxMin scale transformation should occur here - then the scaled results can be loaded into the other dataframe 


    # index and insert column names to normalized data
    result_scaled.index = np.arange(1, len(result) + 1)
    result_scaled.insert(0, 'Percentile', result['Percentile'])
    result_scaled = result_scaled.sort_values("Percentile")

# 2/15/2023 - insert chunk for batch genebody coverage line plot ------------------------------------

    # transpose the columns so they can be group as an x-axis and y-axis   
    result_scaled_noindex = result_scaled.reset_index(drop=True)
    result_scaled_noindex.index = np.arange(1, len(result_scaled_noindex) + 1)
    result_scaled_transposed = result_scaled_noindex.transpose()
    result_scaled_transposed = result_scaled_transposed.reset_index()

# Use the rename() method to make the first row the column headers:
    result_scaled_transposed = result_scaled_transposed.rename(columns=result_scaled_transposed.iloc[0])
    result_scaled_transposed = result_scaled_transposed.drop(result_scaled_transposed.index[0])

# make everything in the dataframe numeric type
    result_scaled_transposed = result_scaled_transposed.apply(pd.to_numeric)

# convert to long (tidy) form
    dfm = result_scaled_transposed.melt('Percentile', var_name='Samples', value_name='Normalized Counts')

    results_tidydata = result_scaled
    results_tidydata.index =  results_tidydata.pop("Percentile")
    results_tidydata = results_tidydata.T

    results_tidydata.to_csv(str(_out_dir)+"/"+str(_project_name)+"_Genebody_Coverage_summary.csv") # the output is the specified out-directory and desired filename
#   remove the first column

if __name__ == "__main__":
### Run Genebody_Coverage_summary_v2.py() to take input using commandline arguments (USER Input)   
    
    parser = argparse.ArgumentParser(description="Genebody Summary script for tab-demlimited un-normalized counts in a sample batch directory -- Argument Parser")

    parser.add_argument("-txtdir", "--txt-dir", type=os.path.abspath, required=True,
                        help="[REQUIRED] Provide the path to the USER input dirctory of tab-delimited txt files for genebodycoverage.\n -txtdir [~/absolute/path/to/tab-delimited-txtdata/directory],\t--txt-dir [~/absolute/path/to/tab-delimited-txtdata/directory]\n")

    parser.add_argument("-out", "--output-directory", required=True, type=os.path.abspath,
                        help="[REQUIRED] Provide the desire output path for the directory containing the plot.\n -out [~/OUTPUT-FILE-PATH/],\t--output-directory [~/OUTPUT-FILE-PATH/]\n")

    parser.add_argument("-name", "--project-name", required=True, 
                        help="[REQUIRED] Name your project for whiich is table is being generated for.")

    args = parser.parse_args()

    _input_txt_dir     = args.txt_dir
    _out_dir         = args.output_directory
    _project_name    = args.project_name
    
    geneBody_summaryTable(_input_txt_dir, _out_dir, _project_name)
