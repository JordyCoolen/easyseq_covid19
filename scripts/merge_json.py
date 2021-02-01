#!/usr/bin/env python

import simplejson as json
import argparse
import os

INFO = "Merge JSON files, where json1 will be used as final final result"

__version__ = 0.1

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json1", type=str, required=True,
                        help="full path to first JSON file")
    parser.add_argument("--json2", type=str, required=True,
                        help="full path to second JSON file")
    parser.add_argument("--key", type=str, required=True,
                        help="Name of the key to add")

    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    args = parser.parse_args()

    return args

def run(args):
    print(f"{args.json1} and {args.json2}")

    # open both the json files and read
    if os.path.exists(args.json1):
        print(f"{args.json1} exists")
        with open(args.json1) as f1:
            data1 = json.load(f1)

    if os.path.exists(args.json2):
        print(f"{args.json2} exists")
        with open(args.json2) as f2:
            data2 = json.load(f2)

    print('test')

    # add one dictionary/json to other json
    data1[args.key] = data2

    # os.path.join(os.path.realpath(args.json1),
    out = open(args.json1, 'w')
    json.dump(data1, out)


if __name__ == "__main__":
    # load arguments
    args = parse_args()
    run(args)
    print("Finished")