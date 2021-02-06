#!/usr/bin/env python

######
INFO = "Convert results to PDF report"
__version__ = 0.3
######

"""
Title:          final_report.py
Author:         J.P.M. Coolen
Date:           06-02-2021 (dd-mm-yyyy)
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

    logo2 = os.path.join(localdir, "report/logo2.png")
    logo2 = logo2.replace(' ','%5')

    lineage_df = pd.read_csv(args.lineage)

    # obtain annotation file and stats
    variant_stats_df = pd.read_csv(args.annotation, sep='\t', engine='python', comment='##')

    # filter only annotation for better overview
    annotation_df = variant_stats_df[['Sample','Position','Var Type','HGVS','Shorthand']]

    # fill html
    template_vars = {
        # pretty things
        "logo": logo,
        "logo2": logo2,
        "version": __version__,

        # general info
        "sampleName": args.sampleName,

        # lineage
        "lineage": lineage_df.to_html(index=False, header=True),

        # annotation
        "annotation": annotation_df.to_html(index=False, header=True),

        # variant stats
        "variant_stats": variant_stats_df.to_html(index=False, header=True),

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
