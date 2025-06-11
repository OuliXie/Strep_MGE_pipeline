#!/usr/bin/env python

import os
import argparse
import shutil
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# Post-processing script for taking a user-defined list of insertion sites and re-orders MGE segments according
# to the user defined list.
# The output of this can then be used cluster MGEs with insertion site data
# Will generate new fasta headers in format <genome>;<insertion site>;<contig location>;<mge class>


def rename_fasta(source, genome, insertion_alias, mge_class):
    # read in fasta - expect only a single fasta sequence
    # will give an exception if give it a multifasta
    record = SeqIO.read(source, "fasta")
    record.id = genome + ";" + insertion_alias + ";" + record.id.split(";")[-1] + ";" + mge_class
    record.description = ""

    return record


def get_insertions(insertions, summary_mge, fasta, output):
    mge_types = ["Phage", "Phage_like", "ICE", "ME"]

    # clean summary_mge dataframe
    summary_mge.replace("IS;", "", regex=True, inplace=True)
    summary_mge.replace(";Hotspot", "", regex=True, inplace=True)
    summary_mge.replace("Unclassified", np.nan, inplace=True)
    summary_mge.replace("IS", np.nan, inplace=True)
    summary_mge.replace("Non_MGE", np.nan, inplace=True)
    summary_mge.replace("Degraded", np.nan, inplace=True)
    # find columns with Phage, Phage_like, ICE or ME
    mge_cols = []
    for col in summary_mge.columns:
        if summary_mge[col].isin(mge_types).sum() > 0:
            mge_cols.append(col)

    # Loop through user insertion sites
    for i in range(len(insertions)):
        pan_insertion = insertions.iloc[i, 0]
        # Check if insertion is a coreless segment
        if "coreless" in pan_insertion:
            # Remove from list of insertion sites as now analysed
            if pan_insertion in mge_cols:
                mge_cols.remove(pan_insertion)
                # Check output folder exists
                coreless_out = os.path.join(output, pan_insertion)
                if not os.path.exists(coreless_out):
                    os.mkdir(coreless_out)
                # Loop through genomes which are positive and copy into output folder
                for genome in summary_mge.loc[summary_mge[pan_insertion].isin(mge_types), "Gff"].tolist():
                    source = os.path.join(fasta, pan_insertion, genome + "_" + pan_insertion + ".fa")
                    destination = os.path.join(coreless_out, genome + "_" + pan_insertion + ".fa")
                    mge_class = summary_mge.loc[summary_mge["Gff"] == genome, pan_insertion].item()
                    # Rename fasta
                    segment = rename_fasta(source, genome, pan_insertion, mge_class)
                    # Write renamed fasta
                    SeqIO.write(segment, destination, "fasta")
        else:
            # Check if name has been designated for this site
            insertion_alias = insertions.iloc[i, 1]
            if insertion_alias not in ["Rearrangement", "rearrangement", "Novel", "novel", ""]:
                insertion_out = os.path.join(output, "i_" + insertion_alias)
                insertion_type = "specified"
            else:
                insertion_alias = pan_insertion
                insertion_type = "new"
                insertion_out = os.path.join(output, insertion_alias)

            # Check if insertion site folder exists
            if not os.path.exists(insertion_out):
                os.mkdir(insertion_out)
            else:
                pass

            # Generate list of possible insertions grouped with user insertion sites
            core_genes = pan_insertion.split("-")
            # Check insertion sites specified correctly
            if len(core_genes) != 2:
                print("Error with insertion site " + pan_insertion +
                      ". Only insertion sites should be separated by \"-\".")
            else:
                pass

            if insertion_type == "specified":
                insertion_list = [pan_insertion, core_genes[1] + "-" + core_genes[0],
                                  "Sequence_break-" + core_genes[0], core_genes[0] + "-Sequence_break",
                                  "Sequence_break-" + core_genes[1], core_genes[1] + "-Sequence_break"]
            else:
                insertion_list = [pan_insertion]

            # Loop through each possible insertion
            for alt_insertion in insertion_list:
                # Check if exact match in mge_cols
                if alt_insertion in mge_cols:
                    # Remove from mge_cols as now analysed
                    mge_cols.remove(alt_insertion)
                    # Loop through genomes with MGE at this insertion
                    for genome in summary_mge.loc[summary_mge[alt_insertion].isin(mge_types), "Gff"].tolist():
                        source = os.path.join(fasta, alt_insertion, genome + "_" + alt_insertion + ".fa")
                        mge_class = summary_mge.loc[summary_mge["Gff"] == genome, alt_insertion].item()
                        # Read in fasta and rename header
                        segment = rename_fasta(source, genome, insertion_alias, mge_class)
                        # Write renamed fasta
                        SeqIO.write(segment, os.path.join(
                            insertion_out, genome + "_" + os.path.basename(insertion_out) + ".fa"), "fasta")

    # Check if any remaining insertion sites not specified
    if len(mge_cols) > 0:
        for j in mge_cols:
            insertion_out = os.path.join(output, j)
            # Check if insertion site folder exists
            if not os.path.exists(insertion_out):
                os.mkdir(insertion_out)
            for genome in summary_mge.loc[summary_mge[j].isin(mge_types), "Gff"].tolist():
                source = os.path.join(fasta, j, genome + "_" + j + ".fa")
                destination = os.path.join(insertion_out, genome + "_" + j + ".fa")
                mge_class = summary_mge.loc[summary_mge["Gff"] == genome, j].item()
                # Rename fasta
                segment = rename_fasta(source, genome, j, mge_class)
                # Write renamed fasta
                SeqIO.write(segment, destination, "fasta")

    return


def main():

    description = "Reorders MGE segments based on user defined list of insertion sites"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--insertion",
        dest="insertion",
        help="csv of user defined insertion sites - col 1 pangenome names, col 2 alias for output",
        required=True,
        type=os.path.abspath
    )

    io_opts.add_argument(
        "-m",
        "--mge",
        dest="mge",
        help="summary_mge_classes.csv from MGE pipeline",
        required=True,
        type=os.path.abspath
    )

    io_opts.add_argument(
        "-f",
        "--fasta",
        dest="fasta",
        help="path to folder of mge segment fastas",
        required=True,
        type=os.path.abspath
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output",
        help="output directory",
        required=False,
        type=os.path.abspath
    )

    args = parser.parse_args()

    # read in summary_mge_classes.csv
    summary_mge = pd.read_csv(args.mge)

    # read in user defined list of insertion sites
    insertions = pd.read_csv(args.insertion, header=None, keep_default_na=False)
    # check insertion csv is in right format
    if len(insertions.columns) != 2:
        print("Error: input insertion site list in wrong format")
        exit(1)

    # check output exists
    if args.output is None:
        output = "ordered_segments"
        if os.path.exists(output):
            print("Default output directory ordered_segments already exists")
            exit(1)
        else:
            os.mkdir("ordered_segments")
    else:
        output = args.output
        if os.path.exists(output):
            print("Specified output directory " + output + " already exists")

    get_insertions(insertions, summary_mge, args.fasta, output)


if __name__ == "__main__":
    main()
