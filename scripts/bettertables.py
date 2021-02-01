#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Title:          bettertables
Company:        Nimagen
Author:         Rick Adriaan Lammerts
Date:           07-01-2021 (dd-mm-yyyy)
Python Version: 3.7.x
Description:    makes better tables from the raw output of bcftools using vcfpy
"""

import argparse
from sys import stdin, stdout, stderr, exit
#from pysam import VariantFile  # pysam was unwieldy to use, switched over to vcfpy
import vcfpy


# this script can use some polish, objects to store the data, the ability to define an output and input file
# it will crash if there is more than one call on a position!
# this was built in because there were no examples available with multiple calls when writing this script


# since the new headers are in a dictionary, the order was defined in a list
# allows the program to only output the headers when merging the files
table_order = [
    "Sample",
    "Position",
    "Var Type",
    "Reference",
    "Alternative",
    "Quality",
    "Total Read Depth",
    "Read Depth Call",
    "Ref Forw",
    "Ref Rev",
    "Alt Forw",
    "Alt Rev"]


def main():
    # get args
    args = get_command_line_arguments()
    
    # print the header
    print("\t".join(table_order))

    # the header can be printed out without running the rest of the program 
    # this was to make joining tables easier while still giving them a header
    if not args.header:
        for vcfdata in vcf2table(args.input, args.sample):
            tabular_out(vcfdata)


def vcf2table(vcfstream, sample_name):
    """turns a vcf file into a dictionary with the new headers, yields every line
    the headers in the dictionary should correspond with the table_order list"""
    
    with vcfpy.Reader.from_stream(vcfstream) as vcffile:
        for rec in vcffile:
            if len(rec.ALT) > 1:  # TODO make compatible with multiple calls per position
                stderr.write("WARNING, more than one ALT!")
                exit(1)
            if sample_name == "CHROM":
                sample_name = rec.CHROM
            yield {
                "Sample":           sample_name,
                "Position":         rec.POS,
                "Var Type":         rec.ALT[0].type,
                "Reference":        rec.REF,
                "Alternative":      rec.ALT[0].value,
                "Quality":          rec.QUAL,
                "Total Read Depth": rec.INFO["DP"],
                "Read Depth Call":  rec.calls[0].data.get("DP"),
                "Ref Forw":         rec.calls[0].data.get("ADF")[0],
                "Ref Rev":          rec.calls[0].data.get("ADR")[0],
                "Alt Forw":         rec.calls[0].data.get("ADF")[1],
                "Alt Rev":          rec.calls[0].data.get("ADR")[1]
            }


def tabular_out(vcfdata, order=table_order):
    """generate output to stdout"""
    print("\t".join([str(vcfdata[item]) for item in order]))


def get_command_line_arguments():
    """Handles the argument in and output on the command line, returns the
    arguments given by the user"""
    argp = argparse.ArgumentParser(description="Reformats vcffiles, only takes plain text vcf files from stdin, outputs table on stdout")
    argp.add_argument("input", default=stdin, nargs="?", type=argparse.FileType('r'),
                      help="""vcf file to parse""")
    argp.add_argument("--sample", type=str, required=False, default="CHROM",
                      help="""sample name for column 1, default=CHROM""")
    argp.add_argument("--header", action="store_true",
                      help="""output header and exit""")
    #argp.add_argument("output", nargs="?", type=str, default=stdout,
    #                  help="""feature not yet functional""")
    return argp.parse_args()


if  __name__ == "__main__":
    main()