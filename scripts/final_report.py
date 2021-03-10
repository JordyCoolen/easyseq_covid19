#!/usr/bin/env python

######
INFO = "Convert results to PDF report"
__version__ = "0.5.1"
######

"""
Title:          final_report.py
Author:         J.P.M. Coolen
Date:           09-03-2021 (dd-mm-yyyy)
Description:    Convert results to PDF report
"""

import os
import argparse
import pandas as pd

def parse_args():
    """
        Argument parser
    """
    parser = argparse.ArgumentParser(description=INFO, \
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sampleName", type=str, required=True,
                        help="name of sequence sample"),
    parser.add_argument("--lineage", type=str, required=True,
                        help="location to lineage file"),
    parser.add_argument("--annotation", type=str, required=False,
                        help="location to annotation file"),
    parser.add_argument("--params", type=str, required=False,
                        help="location to parameters.txt file"),
    parser.add_argument("--stats", type=str, required=False,
                        help="location to stats.txt file"),
    parser.add_argument("--HVdel", type=str, required=False,
                        help="location to kma output file for HV69-70"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def fill_html(args):
    '''
        Code to fill in the placeholders in the html
        and generate a html and pdf

        :params JSON: JSON object containing all the results
        :params outputDir: directory to store the results

        :out pdf: pdf report of the results
        :out html: html report of the results
    '''

    import matplotlib
    matplotlib.use('Agg')
    from weasyprint import HTML
    from jinja2 import Environment, FileSystemLoader

    print('Start Filling')

    localdir = os.path.dirname(os.path.realpath(__file__))

    # create and render html file with tables
    env = Environment(loader=FileSystemLoader(localdir))
    template = env.get_template('report/final_report_template.html')

    # location of logo
    logo = os.path.join(localdir, "report/logo.png")
    logo = logo.replace(' ','%20')

    # load parameters.txt file
    params_df = pd.read_csv(args.params, sep='\t')

    # load parameters.txt file
    stats_df = pd.read_csv(args.stats, sep='\t')
    stats_df = stats_df[['len','N','coverage']]

    # load file with lineage output
    lineage_df = pd.read_csv(args.lineage)

    # obtain annotation file and stats
    variant_stats_df = pd.read_csv(args.annotation, sep='\t', engine='python', comment='##')

    try:
        # calculate ALT freq
        variant_stats_df["ALT freq"] = (variant_stats_df["Alt Forw"]+variant_stats_df["Alt Rev"]) /\
                                       (variant_stats_df["Ref Forw"]+variant_stats_df["Ref Rev"] +
                                       variant_stats_df["Alt Forw"]+variant_stats_df["Alt Rev"]) * 100

        variant_stats_df = variant_stats_df.round({"ALT freq": 1})
    except KeyError:
        variant_stats_df = pd.DataFrame({'NA': []})

    # filter only annotation for better overview
    try:
        annotation_df = variant_stats_df[['Position','Var Type','Read Depth Call',
                                          'ALT freq','HGVS','Shorthand']]
    except KeyError:
        annotation_df = pd.DataFrame({'NA': []})

    # fill html
    template_vars = {
        # pretty things
        "logo": logo,
        "version": __version__,

        # general info
        "sampleName": args.sampleName,

        # lineage
        "lineage": lineage_df.to_html(index=False, header=True),

        # genome stats
        "stats": stats_df.to_html(index=False, header=True),

        # annotation
        "annotation": annotation_df.to_html(index=False, header=True),

        # variant stats
        "variant_stats": variant_stats_df.to_html(index=False, header=True),

        # parameters
        "parameters": params_df.to_html(index=False, header=True),

    }

    # output pdf
    outfile = os.path.join(args.outputDir, '{}.pdf'.format(args.sampleName))

    # render html and write
    html_out = template.render(template_vars)
    with open(os.path.join(args.outputDir, '{}.html'.format(args.sampleName)), 'w') as html_file:
        html_file.write(html_out)

    # save html as pdf to disc
    HTML(string=html_out, base_url=__file__).write_pdf(outfile,
                                                       stylesheets=[os.path.join(localdir, 'report/style.css')])
    
if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    fill_html(args)
    print("Finished")
