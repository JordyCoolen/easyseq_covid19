#!/usr/bin/env python

######
INFO = "Code to extract fastq header information from fastq file"
__version__ = 0.10

######

import argparse
import os
import storage.storage as storage
import logging

logging = logging.getLogger('test')
import pyfastx


def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--inputFile", type=str, required=True,
                        help="input full path of input report file")
    parser.add_argument("--JSON", type=str, required=False,
                        help="full path to JSON file", default=None)
    parser.add_argument("--name", type=str, required=False,
                        help="name to store in JSON as tool name", default='fastq')
    parser.add_argument("--sample", type=str, required=True,
                        help="name of sequence sample"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args


def run(args, name):
    count = 0
    # read fastq line and extract parts
    for read in pyfastx.Fastq(args.inputFile):
        count += 1
        split = read.description.split(':')
        runid = f"{split[0].replace('@', '')}_0{split[1]}_A{split[2]}"
        barcode = split[-1]
        if count == 1:
            break

    # list of tuples of input for JSON
    results = [("runID", runid),
               ("barcode", barcode)]

    # if JSON is present use exiting, else create new unique name
    JSON = storage.JSON()
    if not args.JSON:
        JSON.name(args.sample)
    JSON.open(args.JSON)
    JSON.add_results(args.name, results)
    JSON.pretty_print()
    JSON.write(args.outputDir)

    logging.info(results)


if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    logging.info(args)
    name = args.inputFile.replace('.report', '')
    name = name.replace('.', '_')
    run(args, name)
    logging.info('Finished importing information')
