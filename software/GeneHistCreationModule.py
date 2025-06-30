#!/usr/bin/env python3

"""
Counts Matrix to Gene Histogram

This script reads a tab-delimited raw counts matrix (genes x samples),
transforms counts to log2 CPM (Counts Per Million), and bins the counts
into specified intervals. It outputs a table showing how many genes fall
within each bin for each sample.

Usage:
    python counts_to_genehist.py --rawcounts input_counts.tsv --outputfile genehist.tsv
"""

import argparse
import pandas as pd
import numpy as np

def CountsMatrixToGeneHist(df, binsize=0.25, maxdepth=18.5):
    """
    Converts a raw counts matrix to a binned gene histogram.

    Parameters:
    - df: pandas DataFrame of counts (genes x samples), with optional non-numeric columns in first two columns.
    - binsize: size of each bin interval.
    - maxdepth: maximum bin cutoff.

    Returns:
    - pandas DataFrame with bin ranges and counts per sample.
    """
    # Remove non-numeric columns (first two columns)
    df = df.iloc[:, 2:]

    # Convert to numeric (coerce errors to NaN)
    df = df.apply(pd.to_numeric, errors='coerce')

    # CPM transformation
    read_depth = df.sum(axis=0) / 1e6
    df = df.div(read_depth, axis=1)

    # Log2 transform
    df = np.log2(df + 1)

    # Create bins
    bins = np.arange(0, maxdepth, binsize)

    # Count genes per bin
    bin_counts = []
    for bin in bins:
        top_bin = bin + binsize
        count_df = (df > bin) & (df <= top_bin)
        count_df_sum = count_df.sum(axis=0)
        bin_counts.append(count_df_sum)

    # Assemble output DataFrame
    bin_counts_df = pd.DataFrame(bin_counts)
    topbins = bins + binsize
    xaxis_names = [f"({bin},{topbin}]" for bin, topbin in zip(bins, topbins)]
    bin_counts_df.insert(loc=0, column='Bins', value=xaxis_names)

    return bin_counts_df


def main():
    parser = argparse.ArgumentParser(
        description="Convert a raw counts matrix to a binned gene histogram table."
    )
    parser.add_argument(
        "--rawcounts",
        required=True,
        help="Path to raw counts input file (tab-delimited)."
    )
    parser.add_argument(
        "--outputfile",
        required=True,
        help="Path to output TSV file with bin counts."
    )

    args = parser.parse_args()

    # Read counts matrix
    df = pd.read_csv(args.rawcounts, sep="\t")

    # Process
    bin_counts_df = CountsMatrixToGeneHist(df)

    # Write output
    bin_counts_df.to_csv(args.outputfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
