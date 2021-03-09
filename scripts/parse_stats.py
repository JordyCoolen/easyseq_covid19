#!/usr/bin/env python

######
INFO = "Small parser to calculate stats"
######

"""
Title:          parse_stats.py
Author:         J.P.M. Coolen
Date:           09-03-2021 (dd-mm-yyyy)
Description:    Small parser to calculate stats
"""

import argparse
import pandas as pd

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stats", type=str, required=True,
                        help="location to stats file"),

    # parse all arguments
    args = parser.parse_args()

    return args

def main(args):

    # obtain annotation file and stats
    df = pd.read_csv(args.stats, sep='\t', engine='python',  index_col=0)

    # remove total row
    try:
        df = df.drop(index='total')
    except KeyError:
        pass

    # calculate coverage
    df['coverage'] = round(((29903 - df['N'])/29903) * 100, 2)

    # overwrite to stats file including coverage
    df.to_csv(args.stats, sep='\t')

if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    main(args)