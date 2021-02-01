#!/usr/bin/env python

######
import sys

INFO = "Convert results to PDF report"
__version__ = 0.2
######

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
    parser.add_argument("--variants", type=str, required=True,
                        help="location to variants file"),
    parser.add_argument("--genes", type=str, required=True,
                        help="location to genes mutation file"),
    parser.add_argument("-o", "--outputDir", type=str, required=False,
                        help="full path of output folder", default=os.path.abspath("./"))
    parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s {version}".format(version=__version__))

    # parse all arguments
    args = parser.parse_args()

    return args

def check_contents(data):
    """
        code to check if expected keys are present in JSON
        :param JSON: JSON object containing all the results
        :return:
    """

    exp_keys = {'fastq': [],
                'centrifuge': [],
                'contaminants': [],
                'KMA': [],
                'WGS_typing': [],
                'variants': ['resistance',
                             'resistance_gene',
                             'all_variants',
                             'unknown',
                             'all_drugs'],
                'BAM_QC': ['avg_dp',
                           'low_cov_frac',
                           'dp_thres',
                           'cov_thres',
                           'qc_passed']}
    try:
        for k in exp_keys:
            data[k]
            for s in exp_keys[k]:
                data[k][s]
    except KeyError as e:
        key = e.args[0]
        print('key {} not present in json, skip report generator'.format(key))
        exit(0)


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

    lineage_df = pd.read_csv(args.lineage)
    variants_df = pd.read_csv(args.variants, sep='\t', engine='python', comment='##')
    annotation_df = pd.read_csv(args.annotation, sep='\t', engine='python', comment='##')
    genes_df = pd.read_csv(args.genes, sep='\t', engine='python', comment='# ')

    # fill html
    template_vars = {
        # pretty things
        "logo": logo,
        "version": __version__,

        # general info
        "sampleName": args.sampleName,

        # lineage
        "lineage": lineage_df.to_html(index=False, header=True),

        # variants
        "variants": variants_df.to_html(index=False, header=True),

        # annotation
        "annotation": annotation_df.to_html(index=False, header=True),

        # annotation
        "genes": genes_df.to_html(index=False, header=True),

        # "runID": data['fastq']['runID'],
        # "barcode": data['fastq']['barcode'],
        #
        # # data tables
        # "CENTRIFUGE": CF_df.to_html(index=True, header=False),
        # "CENTRIFUGEDATA": CFd_df.to_html(index=False, header=True),
        # "BAMQC": BAM_QC.to_html(header=False),
        # "CONTAMINANTS": CON_df.to_html(index=False, header=True),
        # "SNPIT": SNPIT.to_html(index=True, header=False),
        # "HSP": HSP.to_html(index=True, header=False),
        # "WGS_IDEN": WGS_IDEN.to_html(index=False, header=True),
        # "SUSCEPTIBILITY": SUS.to_html(index=False),
        # "RESISTANCE": RES_df.to_html(index=False, header=True),
        # "RGENE": RGENE_df.to_html(index=False, header=True),
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

    sys.exit(0)

    # OLD

    # load contents from JSON
    data = JSON.data
    # check if expected keys are present
    check_contents(data)

    # centrifuge
    CF_df = pd.DataFrame.from_dict(data['centrifuge'], orient='index')
    CF_df = CF_df.drop(['metric', 'data', 'barplot'])

    # centrifuge data
    CFd_df = pd.DataFrame.from_dict(data['centrifuge']['data'], orient='index')
    CFd_df = CFd_df.drop(['taxRank'])
    CFd_df = CFd_df.transpose()

    #BAM QC
    BAM_QC = pd.DataFrame.from_dict(data['BAM_QC'], orient='index')

    # contaminants
    CON_df = pd.DataFrame.from_dict(data['contaminants']['data'], orient='index')
    CON_df = CON_df.drop(['taxRank'])
    CON_df = CON_df.transpose()

    # Identification methods
    if "SNP-IT" in data.keys():
        SNPIT = pd.DataFrame.from_dict(data['SNP-IT'], orient='index')
        SNPIT = SNPIT.drop(['sample_name'])
    else:
        SNPIT = pd.DataFrame(['NO DATA'])

    # HSP65 identification
    HSP = pd.DataFrame.from_dict(data['KMA'], orient='index')
    HSP = HSP.drop(['data'])

    # WGS typing
    WGS_IDEN = pd.DataFrame.from_dict(data['WGS_typing']['data'], orient='index')
    WGS_IDEN = WGS_IDEN.drop(['taxRank'])
    WGS_IDEN = WGS_IDEN.drop(['taxID'])
    WGS_IDEN = WGS_IDEN.transpose()

    # get detected resistances
    RES = data['variants']['resistance']

    # create pandas dataframe of resistance mutation
    RES_df = get_resistance_mutations(RES)

    # obtain pandas dataframe of target gene mutations
    RGENE = data['variants']['resistance_gene']
    RGENE_df = get_resistance_gene(RGENE)

    # create susceptibility table
    SUS = susceptibility(RES_df, data)
    
if __name__ == "__main__":
    # load arguments set global workdir
    args = parse_args()
    fill_html(args)
    print("Finished")
