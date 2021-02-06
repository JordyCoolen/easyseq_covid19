#!/usr/bin/env python

"""
Title:          HV69-70
Author:         J.P.M. Coolen
Date:           06-02-2021 (dd-mm-yyyy)
Python Version: 3.7.x
Description:    Parses the KMA HV 69-70 deletion output
"""

import pandas as pd
import argparse
import sys
import os

def get_command_line_arguments():
    """Handles the input KMA results file and output on the command line, returns the
    arguments given by the user"""
    argp = argparse.ArgumentParser(description="Code to parse KMA output for the HV 69-70 deletion")
    argp.add_argument("-i", "--input", type=str, required=True,
                      help="""KMA output file to parse""")
    argp.add_argument("-p", "--path", type=str, required=True,
                      help="""full path to the vcf files""")
    return argp.parse_args()

def main():
    """Steps of the scripts in main"""
    # get args
    args = get_command_line_arguments()
    # parse KMA results file
    result = read_kma_res(args.input)

    if result == 'Wildtype (NC_045512.2)':
        with open(os.path.join(args.path,'WT.vcf'), "r") as output:
            sys.stdout.write(output.read())
    elif result == 'HV 69-70 deletion':
        with open(os.path.join(args.path,'HV69-70.vcf'), "r") as output:
            sys.stdout.write(output.read())

def read_kma_res(input):
    """Read the KMA results of the HV 69-70 detection"""
    # read kma HV 69-70 results
    HVdel_df = pd.read_csv(input, sep='\t', engine='python', comment='##')
    HVdel_df = HVdel_df.sort_values(by=['Score'], ascending=False)
    HVdel_df = HVdel_df.reset_index(drop=True)
    HVdel = HVdel_df.iloc[0][0]
    return HVdel

if  __name__ == "__main__":
    """main"""
    main()