#!/usr/bin/env python

import pandas as pd
import argparse
import os
from Bio import SeqIO
from datetime import datetime

INFO = "Script to convert fasta to GISAID id according to dataTable"
__version__ = '0.1'

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--dataTable", type=str, required=True,
                        help="full path to dataTable")
    parser.add_argument("--overviewFile", type=str, required=True,
                        help="full path to dataTable")
    parser.add_argument("-f", "--fastaFile", type=str, required=True,
                        help="input full path of multi-fasta file")
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

# read datafile
def read_dataTable(file):
    df = pd.read_excel(file, index_col=0, dtype={'MNR': str})
    return(df)

# read fasta files
def match_samples(fastaFile, outputDir, df):

    matches = []
    renamed = []
    log = []

    with open(fastaFile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for index, sample in enumerate(df['MNR']):
                #change name to name before "_"
                record.id = record.id.split('_')[0]
                if record.id == sample:
                    print(f"match: {record.id} --> {df['Stamnaam'].iloc[index]}")
                    log.append(f"{record.id}\t{df['Stamnaam'].iloc[index]}")
                    record.id = df['Stamnaam'].iloc[index]
                    record.description = ''

                    # check if present
                    if record.id not in matches:
                        matches.append(record.id)
                        renamed.append(record)
                    else:
                        for i, entry in enumerate(renamed):
                            if entry.id == record.id:
                                if len(record.seq) > len(entry.seq):
                                    print(f'better hit: {record.id}')
                                    log.append(f"better hit:{record.id}")
                                    renamed[i] = record

    time = datetime.now().strftime("%Y%m%d")
    SeqIO.write(renamed, os.path.join(outputDir, f"{time}_GISAID.fasta"), "fasta")
    with open(os.path.join(outputDir, f"{time}_GISAID.log"), 'w') as logout:
        for l in log:
            logout.write(f"{l}\n")

    print(f"Length: {len(matches)}")
    print(f"Length: {len(set(matches))}")
    return(matches)

# read fasta files
def match_samples2(fastaFile, df):
    match = []
    #matches = []
    #renamed = []

    for index, sample in enumerate(df['MNR']):
        with open(fastaFile) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id == str(sample):
                    print(record.id)
                    match.append(record)

# convert ID's
def get_GISAID(matches, fastaFile, df):
    df = df[df['MNR'].isin(matches)]

    with open(fastaFile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for index, match in enumerate(matches):
                if record.id == match:
                    print(index)
                    print(match)
                    print(df["Stamnaam"].iloc[index])
                    break

    # drop MNR column

# export fasta files with correct ID's

# filter for >=90% coverage
def filter_90(file):
    # read excel overview file
    df = pd.read_excel(file, dtype={'MNR': str})
    # filter >= 90
    df = df.loc[df['coverage'] >= 90]
    #extract Monster number
    MNR = df['MNR'].tolist()
    return(MNR)

def filter_dataTable(df, MNR):
    # extract only samples in MNR
    df = df.loc[df['MNR'].isin(MNR)]
    return(df)

def run(args):
    print('run')
    df = read_dataTable(args.dataTable)
    print(df)
    MNR = filter_90(args.overviewFile)
    print(MNR)
    df = filter_dataTable(df, MNR)
    print(df)
    matches = match_samples(args.fastaFile, args.outputDir, df)

if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    run(args)