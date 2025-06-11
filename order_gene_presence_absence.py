#!/usr/bin/env python

import os
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm

# Script for sorting annotated_gene_presence_absence
# Uses associated gff to fix issues with re found genes and out of order locus tags from PGAP


def parse_gff(file):
    # Adapted from post_run_gff_output.py
    # Read in and split off FASTA portion
    raw_file = open(file, 'r')
    # Read in and split off FASTA portion
    lines = raw_file.read().replace(',', '')
    split = lines.split('##FASTA')[0]
    # Remove headers
    body = []
    for line in split.splitlines():
        if "##" not in line:
            body.append(line)
    # Parse locus tag only in order and save associated gff interval
    parsed_gff = {}
    chrom_list = {}
    for gff_line in body:
        initial_split = gff_line.split("\t")
        attribute_split = " ".join(initial_split[8:]).split(";")
        chromosome = initial_split[0]
        start = initial_split[3]
        end = initial_split[4]
        # Get order in chromosome
        if chromosome in chrom_list.keys():
            chrom_list[chromosome] += 1
            order = chrom_list[chromosome]
        else:
            chrom_list[chromosome] = 1
            order = 1
        locus_tag = ""
        for attribute in attribute_split:
            if "locus_tag" in attribute:
                locus_tag = attribute.split("=")[1]
        parsed_gff[locus_tag] = [chromosome, start, end, order]

    return parsed_gff


def order_extract(genome, gpa, gff_list, out_dir):
    # Read in gff
    gff_path = [x for x in gff_list if os.path.basename(x).startswith(genome + ".") or
                os.path.basename(x).startswith(genome + "_")]
    # Check only one gff file found
    if len(gff_path) == 1:
        gff_path = gff_path[0]
    elif len(gff_path) == 0:
        print("Error: No gff found for " + str(genome))
        exit(1)
    else:
        print("Error: " + str(len(gff_path)) + " files found for " + "genome")
        print(gff_path)
        exit(1)
    gff_locus = parse_gff(gff_path)
    # Select columns from annotated_gene_presence_absence.csv
    gpa = gpa[["Gene", "recombinase", "T4SS", "phage", "Non-unique Gene name", "Annotation", "No. isolates",
                       "No. sequences", genome]]
    # Temp column containing only the first locus_tag for CDS which were found to belong to the same gene
    unique_tag = [x.split(";")[0] for x in gpa[genome].to_list()]
    gpa = gpa.assign(locus_unique=unique_tag)
    # Check which locus tags have been merged in gff
    gff_unique = [x for x in gff_locus.keys() if x in unique_tag]
    # Subset annotated_gene_presence_absence.csv
    gpa = gpa[gpa.loc[:, "locus_unique"].isin(gff_unique)]
    gpa = gpa.set_index("locus_unique")
    # Reorder based on gff order
    gpa = gpa.reindex(gff_unique)

    # Add chromosome and gene start and end positions (gff format)
    cds_list = [x.split(";") for x in gpa[genome].to_list()]
    start = []
    end = []
    chrom = []
    # Iterate through each gene
    for i in cds_list:
        interval_start = []
        interval_end = []
        chrom_list = []
        order_list = []
        for j in i:
            chrom_list.append(gff_locus[j][0])
            interval_start.append(gff_locus[j][1])
            interval_end.append(gff_locus[j][2])
            order_list.append(gff_locus[j][3])
        # Check the CDS are on the same chromosome and if so, are adjacent
        if (len(set(chrom_list)) == 1) and (max(order_list) - min(order_list) < len(chrom_list)):
            start.append(min(interval_start))
            end.append(max(interval_end))
            chrom.append(chrom_list[0])
        # Some genes on the ends of contigs are re-found and merged across chromosomes
        # Because will be ordered by the first CDS, take the interval and contig of that CDS
        # Should not affect clustering in the end as will only be missing small amount of genetic information
        else:
            start.append(interval_start[0])
            end.append(interval_end[0])
            chrom.append(chrom_list[0])
    # Add to gpa
    gpa = gpa.assign(chromosome=chrom, start=start, end=end)

    # Remove phage annotations if gene is also a recombinase
    gpa["phage"] = np.where(gpa["recombinase"] != "", "", gpa["phage"])

    gpa.to_csv(out_dir + "/" + genome + "_gene_presence_absence.csv", index=False)

    return


def main():
    description = "Orders annotated_gene_presence_absence.csv by each genome"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_file",
        help="annotated_gene_presence_absence.csv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-g",
        "--gffs",
        dest="gffs",
        help="Path to post-Panaroo corrected gff folder",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        help="Output directory",
        required=False,
        type=os.path.abspath,
        default="ordered_gene_presence_absence"
    )

    args = parser.parse_args()

    # Make dedicated directory for outputs
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # Read in annotated_gene_presence_absence
    gpa = pd.read_csv(args.input_file, low_memory=False).fillna("")

    # Make copy in output folder
    gpa.to_csv(args.output_dir + "/" + "annotated_gene_presence_absence.csv", index=False)

    # Get list of genomes
    exclude = ["Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences",
               "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment",
               "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc",
               "Avg group size nuc", "recombinase", "T4SS", "phage"]
    genomes = gpa.loc[:, ~gpa.columns.isin(exclude)].columns

    # Get list of gffs
    gff_list = [os.path.join(args.gffs, x) for x in os.listdir(args.gffs) if x.endswith('.gff')]

    for genome in tqdm(genomes):
        order_extract(genome, gpa, gff_list, args.output_dir)

    # Write sequence_id.txt file
    pd.DataFrame(genomes.to_list()).to_csv("sequence_id.txt", index=False, header=False)


if __name__ == "__main__":
    main()
