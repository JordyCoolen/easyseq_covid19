#!/usr/bin/env python

import pandas as pd
import glob
import argparse
import os
import sys
import datetime

def arg_parser():
    """Handles the argument in and output on the command line, returns the
    arguments given by the user"""
    argp = argparse.ArgumentParser(description="Aggregates all the annotation tables to a single overview")
    argp.add_argument("-i", "--inputDir", default=str,
                      help="""location containing all annot_table.txt files to parse""")
    argp.add_argument("-o", "--outputDir", default=str, nargs="?",
                      help="""output file of the program, default = stdout""")
    return argp.parse_args()

def get_file_paths(inputDir, folder):
    ls = glob.glob(f"{inputDir}/*/{folder}")

    return ls

def check_content(ls):

    # result of paths
    dict = {'lineages': [], 'stats': [], 'mutations': []}

    for path in ls:

        lineage_path = os.path.join(path, 'lineages.txt')
        stats_path = os.path.join(path, 'stats.txt')
        mutation_path = os.path.join(path, 'annotation', 'AGGREGATION_annot_table_SH.txt')

        if os.path.exists(lineage_path):
            dict['lineages'].append(lineage_path)
            pass
        else:
            sys.exit(f"Error: {lineage_path} does not exist")

        if os.path.exists(stats_path):
            dict['stats'].append(stats_path)
            pass
        else:
            sys.exit(f"Error: {stats_path} does not exist")

        if os.path.exists(mutation_path):
            dict['mutations'].append(mutation_path)
            pass
        else:
            sys.exit(f"Error: {mutation_path} does not exist")

    return(dict)

def merge_lineage(dict):
    dataframes = []
    for f in dict['lineages']:
        df = pd.read_csv(f, sep=',', engine='python', comment='##')
        dataframes.append(df)

    appended_data = pd.concat(dataframes)
    print(appended_data.index)

    # remove extra headers
    appended_data = appended_data[appended_data.lineage != 'lineage']
    appended_data = appended_data.drop_duplicates('taxon', keep='first')
    appended_data = appended_data.set_index('taxon')

    print(appended_data)

    return(appended_data)

def merge_stats(dict):
    dataframes = []
    for f in dict['stats']:
        df = pd.read_csv(f, sep='\t', engine='python', comment='##')
        dataframes.append(df)

    appended_data = pd.concat(dataframes)
    print(appended_data)

    # remove extra headers
    appended_data = appended_data[appended_data["#seq"] != '#seq']

    appended_data = appended_data.drop_duplicates('#seq', keep='first')
    appended_data = appended_data.set_index('#seq')

    return(appended_data)

def merge_mutations(dict):
    dataframes = []
    for f in dict['mutations']:
        df = pd.read_csv(f, sep='\t', engine='python', comment='##')
        dataframes.append(df)

    appended_data = pd.concat(dataframes, axis=0)

    appended_data = appended_data.drop_duplicates('Sample', keep='last')
    appended_data = appended_data.set_index('Sample')

    print(appended_data)
    return(appended_data)

def summarize(lineages, stats, mutations):

    return(pd.concat([lineages, stats, mutations], axis=1))

def main():

    date = datetime.datetime.now().strftime('%Y%m%d')
    args = arg_parser()
    ls = get_file_paths(args.inputDir, 'summary')
    dict = check_content(ls)

    lineages_df = merge_lineage(dict)
    stats_df = merge_stats(dict)
    mutations_df = merge_mutations(dict)

    # lineages_df.to_excel(f'{date}_lineage.xlsx')
    # stats_df.to_excel(f'{date}_stats.xlsx')
    # mutations_df.to_excel(f'{date}_mutations.xlsx')

    summary_df = summarize(lineages_df, stats_df, mutations_df)

    # sort the columns on sample
    summary_df.sort_index(axis=1)

    # save aggregation
    summary_df.to_excel(f'{args.outputDir}/{date}_COVID19_summary.xlsx', index_label="MNR")
    summary_df.to_csv(f'{args.outputDir}/{date}_COVID19_summary.csv', sep="|", index_label="MNR")

if __name__ == "__main__":
    # load arguments set global workdir
    # args = parse_args()
    # fill_html(args)
    print("Start")
    main()
    print("Finished")