#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Title:          vcf2table
Company:        Nimagen
Author:         Rick Adriaan Lammerts
Date:           07-01-2021 (dd-mm-yyyy)
Python Version: 3.7.x
Description:    makes more human readable tables from vcf files using vcfpy
"""

import re
import vcfpy
import argparse
from collections import OrderedDict
from sys import stdin, stdout, stderr, exit


# GLOBALS
# error code dict:
ERROR_CODES = {
    "MultipleCall_Err": 1,
    "snpEff_Err": 2
}

MISSENSE_SNPEFF_TAG = "missense_variant"

AA_CODES = {
    "Ala":	"A",
    "Arg":	"R",
    "Asn":	"N",
    "Asp":	"D",
    "Cys":	"C",
    "Glu":	"E",
    "Gln":	"Q",
    "Gly":	"G",
    "His":	"H",
    "Hyp":	"O",
    "Ile":	"I",
    "Leu":	"L",
    "Lys":	"K",
    "Met":	"M",
    "Phe":	"F",
    "Pro":	"P",
    "Glp":	"U",
    "Ser":	"S",
    "Thr":	"T",
    "Trp":	"W",
    "Tyr":	"Y",
    "Val":	"V",
    "*":    "*"
}


def main():
    # get args
    args = get_command_line_arguments()

    #with open(args.input, "r") as infile:
    for n, vcfdata in enumerate(vcf2table(args.input, args.sample, 
                                args.annotations, args.missense_only, args.remove_duplicate_snpeff)):
        if n == 0:
            header_out(vcfdata, args.output)    
        pass
        tabular_out(vcfdata, args.output)


def vcf2table(vcfstream, sample_name, parse_annotations, missense_only, remove_duplicate_snpeff):
    """turns a vcf file into an ordered dictionary with the new headers, yields every line
    the headers in the dictionary should correspond with the table_order list
    Errors out if there is more than one call on every location, this still needs to be changed"""
    HGVS_parse_regx = get_hgvs_parse_patterns()  # for regex efficency
    with vcfpy.Reader.from_stream(vcfstream) as vcffile:
        for rec in vcffile:
            
            if len(rec.ALT) > 1:  # TODO make compatible with multiple calls per position
                stderr.write("WARNING, more than one ALT!\nExiting\n")
                exit(ERROR_CODES["MultipleCall_Err"])
           
            if sample_name == "CHROM":
                sample_name = rec.CHROM
           
            if parse_annotations:
                for annotations in parse_snpeff_annotations(rec, missense_only, remove_duplicate_snpeff, HGVS_parse_regx):
                    vcfdict = get_default_vcf_parse_dict(rec, sample_name)
                    vcfdict["HGVS"], vcfdict["Shorthand"] = annotations
                    yield vcfdict
            else:
                vcfdict = get_default_vcf_parse_dict(rec, sample_name)
                yield vcfdict


def get_default_vcf_parse_dict(rec, sample_name):
    """Gets called by vcf2table for getting all the required vcf fields
    More are added by vcf2table if certain critera are met"""
    return OrderedDict([
                ("Sample",           sample_name),
                ("Position",         rec.POS),
                ("Var Type",         rec.ALT[0].type),
                ("Reference",        rec.REF),
                ("Alternative",      rec.ALT[0].value),
                ("Quality",          rec.QUAL),
                ("Total Read Depth", rec.INFO["DP"]),
                ("Read Depth Call",  rec.calls[0].data.get("DP")),
                ("Ref Forw",         rec.calls[0].data.get("ADF")[0]),
                ("Ref Rev",          rec.calls[0].data.get("ADR")[0]),
                ("Alt Forw",         rec.calls[0].data.get("ADF")[1]),
                ("Alt Rev",          rec.calls[0].data.get("ADR")[1])
            ])


# snpEff header:
# 0        1            2                   3           4         5              6            7                    8       9        10      11                       12                     13                   14         15                   
# 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS
def parse_snpeff_annotations(rec, missense_only, remove_duplicate_snpeff, HGVS_parse_regx):
    """Top level function to parse HGVS annotations, has option to select for only missense variants
    These are checked against a global variable. yields dna annotation if there is no protein annotation.
    Has an option to skip duplicate annotations, as a new dictionary will be made for each annotation
    Needs a precompiled regex dictionary for efficiency
    currently only has shorthands for aa substitutions"""
    try:
        annotation_mem = []
        for annotation in rec.INFO["ANN"]:
            annotation = annotation.split("|")
            if (annotation[1] == MISSENSE_SNPEFF_TAG and missense_only) or not missense_only:
                if annotation[10]:
                    HGVS = ":".join([annotation[3], annotation[10]])
                    Shorthand = hgvs_aasub2shorthand(HGVS, HGVS_parse_regx["aa_sub"])
                else: 
                    # if there is no protein notation yield the DNA one
                    HGVS = ":".join([annotation[3], annotation[9]])
                    Shorthand = "-"
                if (remove_duplicate_snpeff and not HGVS in annotation_mem) or not remove_duplicate_snpeff:
                    annotation_mem.append(HGVS)
                    yield HGVS, Shorthand
    except Exception as snpEff_Error:
        stderr.write(f"snpEff annotation_error!\n{snpEff_Error}\nExiting\n")
        exit(ERROR_CODES["snpEff_Err"])


def get_hgvs_parse_patterns():
    """stores the regex search patterns for parsing HGVS for efficiency.
    Currently only contais the regex for aa substitutions
    Uses the amino acid list so that the pattern is as specific as possible e.g. Ala|Arg e.c.t.
    aasub:
        group 1 = old AA
        group 2 = coordinate
        group 3 = new AA"""
    aa_list = "|".join(AA_CODES.keys())
    aa_list = aa_list.replace(r"*", r"\*")
    return {
        "aa_sub":  re.compile(rf"^\w+:p\.({aa_list})(\d+)({aa_list})$")  # searches for <Gene_Name>:p.<three letter AA list><coordinates><three letter AA list>
    }


def hgvs_aasub2shorthand(hvgs, regx_pattern):
    """Converts the HGVS of a substitution to a shorthand annotation, uses the regex made in get_hgvs_parse_patterns"""
    parsed_hvgs = regx_pattern.match(hvgs)
    if parsed_hvgs:
        # makes the list of 1 letter AA code(regx group 1), coordinate (regx group 2) and 1 letter AA code(regx group 3)
        return "".join([
            AA_CODES[parsed_hvgs[1]],
            parsed_hvgs[2],
            AA_CODES[parsed_hvgs[3]]
        ])
    else:
        return "-"


def header_out(vcfdata, outfile):
    """Print the header from the ordered dictionary"""
    outfile.write("\t".join([str(item) for item in vcfdata.keys()])+"\n")


def tabular_out(vcfdata, outfile):
    """generate output to stdout"""
    outfile.write("\t".join([str(item) for item in vcfdata.values()])+"\n")
    #print("\t".join([str(vcfdata[item]) for item in order]))


def get_command_line_arguments():
    """Handles the argument in and output on the command line, returns the
    arguments given by the user"""
    argp = argparse.ArgumentParser(description="Reformats vcffiles, only takes plain text vcf files from stdin, outputs table on stdout WARNING! will error out when more than one ALT is present!")
    argp.add_argument("input", default=stdin, nargs="?", type=argparse.FileType('r'),
                      help="""vcf file to parse""")
    argp.add_argument("--sample", type=str, required=False, default="CHROM",
                      help="""sample name for column 1, default=CHROM""")
    argp.add_argument("-a", "--annotations", action="store_true",
                      help="""parse snpEff annotations""")
    argp.add_argument("-d", "--remove_duplicate_snpeff", action="store_true",
                      help="""snpeff has a tendency to report duplicates, this is most likely caused by the covid genome having two
                      higly similar/ the same features that were annotated. say GU280_gp01 and GU280_gp01.2""")
    argp.add_argument("-m", "--missense_only", action="store_true",
                      help="""Select only missense mutations in annotations. labled by snpEff as \"missense_variant\"""")
    argp.add_argument("-o", "--output", default=stdout, nargs="?", type=argparse.FileType('w'),
                      help="""output file of the program, default = stdout""")
    return argp.parse_args()


if  __name__ == "__main__":
    main()
