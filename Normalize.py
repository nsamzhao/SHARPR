#!/usr/bin/env python

"""SHARPR Normalize Module

Take RNA and DNA counts from a MPRA experiment and return a normalized log2
ratio for each unique sequence

usage: Normalize.py [-h] [--pcr PCR] [--pcd PCD] [--cutoff CUTOFF]
                    rnalist dnalist outfile

positional arguments:
  rnalist          Specify the name of RNA files separated by ,
  dnalist          Specify the corresponding DNA files separated by ,
  outfile          specify the output file name

optional arguments:
  -h, --help       show this message and exit
  --pcr PCR        (int) pseudo RNA counts added for smoothing
  --pcd PCD        (int) pseudo DNA counts added for smoothing
  --cutoff CUTOFF  (int) threshold of the min number of DNA counts to be
                   retained in the data the default number is 20

"""


__author__ = 'Nanxiang(Samuel) Zhao'
__email__ = 'samzhao@umich.edu'


import numpy as np
import pandas as pd
import argparse
import itertools


# command line input settings
parser = argparse.ArgumentParser(description='Take RNA and DNA counts from a MPRA experiment and ' \
                                             'return a normalized log2 ratio for each unique sequence')
parser.add_argument("rnalist", help='Specify the name of RNA files separated by ,')
parser.add_argument("dnalist", help='Specify the corresponding DNA files separated by ,')
parser.add_argument("outfile", help='specify the output file name')
parser.add_argument("--pcr", help='(int) pseudo RNA counts added for smoothing', default=1, type=int)
parser.add_argument("--pcd", help='(int) pseudo DNA counts added for smoothing', default=1, type=int)
parser.add_argument("--cutoff", help='(int) threshold of the min number of DNA counts to be retained in the data ' \
                                     'the default number is 20', default=20, type=int)

args = parser.parse_args()
rnalist = args.rnalist.split(',')
dnalist = args.dnalist.split(',')


# loop through corresponding rna and dna files together
for rnafiles, dnafiles in itertools.izip_longest(rnalist, dnalist):

    # loading files
    df_rna = pd.read_csv(rnafiles, sep='\t', index_col=0)
    df_dna = pd.read_csv(dnafiles, sep='\t', index_col=0)
    df_dna_copy = df_dna.copy()

    # check pcr and pcd input format
    if not isinstance(args.pcd, int) & isinstance(args.pcr, int):
        print('WARNING: pseudo counts are set as default value 1 since input of pcr and/or pcd is not integer')
    else:
        df_rna += int(args.pcr)
        df_dna += int(args.pcd)

    # computing the normalized log2 ratio
    # the RNA and DNA counts are first divided by their total counts in the experiment
    # then log2(normalized_rna) - log2(normalized_dna) to compute the ratio
    df_rna_sum = df_rna.sum()
    df_dna_sum = df_dna.sum()
    df_rna_normalized = np.log2(df_rna.div(df_rna_sum))
    df_dna_normalized = np.log2(df_dna.div(df_dna_sum))
    df_normalized = pd.DataFrame(df_rna_normalized - df_dna_normalized)

    # filter out counts which are less than the cutoff values in dna file
    df_normalized = df_normalized[df_dna_copy[:] >= args.cutoff]

    # compute the median of the ratio and calculate the number of barcode
    df_normalized['Ratio'] = df_normalized.median(axis=1, skipna=True)
    df_normalized['#barcodes'] =df_normalized.loc[:, df_normalized.columns != 'Ratio'].count(axis=1)

    # output the file
    df_output = pd.concat([df_normalized['Ratio'], df_normalized['#barcodes']], axis=1)
    df_output.to_csv(args.outfile,sep='\t',encoding='utf-8')

