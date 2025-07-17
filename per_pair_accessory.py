#!/usr/bin/env python

import os
import sys
from contextlib import contextmanager
import argparse
import pandas as pd
import numpy as np
import pybedtools
import subprocess


# Script for extracting user specified accessory genes between core-core segments
# If segment contains a recombinase gene, it will also extract the fasta and gff of the segment
# Run order_annotated_gene_presence_absence.sh first in same directory
# order_annotated_gene_presence_absence.sh will also generate a sequence_id.txt file listing all genomes

# If there is read-through to core genes, these will be automatically removed


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
            body.append(line.replace(" ", "_"))
    gff_file = pybedtools.BedTool("\n".join(body), from_string=True)

    return gff_file


def getfasta(bedtools_start, bedtools_end, chrom, sequence, output, fasta_list, gff_list):
    # Make folder for fasta files
    fasta_folder = os.path.join(os.path.dirname(output), "segment_fasta", os.path.basename(output))
    if not os.path.exists(fasta_folder):
        os.mkdir(fasta_folder)
    # Make folder for gff files
    gff_folder = os.path.join(os.path.dirname(output), "segment_gff", os.path.basename(output))
    if not os.path.exists(gff_folder):
        os.mkdir(gff_folder)
    # Set interval for extraction in bed format - note these coordinates are already in bed format
    bed = pybedtools.BedTool(str(chrom) + "\t" + str(bedtools_start) + "\t" + str(bedtools_end) + "\t" +
                             sequence + ";" + os.path.basename(output) + ";" + str(chrom) + ":" +
                             str(bedtools_start) + "-" + str(bedtools_end), from_string=True)
    fasta_file = [x for x in fasta_list if os.path.basename(x).startswith(sequence + ".") or
                  os.path.basename(x).startswith(sequence + "_")]
    # Read in gff
    gff_file = [x for x in gff_list if os.path.basename(x).startswith(sequence + ".") or
                os.path.basename(x).startswith(sequence + "_")]
    # Check only one match
    if len(fasta_file) == 1 and len(gff_file) == 1:
        fasta_file = fasta_file[0]
        # Need to split off FASTA from gff to read into bedtools
        gff_file = parse_gff(gff_file[0])
    else:
        print("Error: " + str(len(fasta_file)) + " matches for " + sequence + ".fa or " + sequence + ".fasta and " +
              str(len(gff_file)) + " matches for gff.")
        exit(1)
    # Check if fasta has been index
    if not os.path.exists(fasta_file + ".fai"):
        faidx_cmd = ["samtools", "faidx", fasta_file]
        subprocess.run(faidx_cmd, check=True)
    # Extract using bedtools getfasta
    fasta = bed.sequence(fi=fasta_file, nameOnly=True)
    gff = pybedtools.BedTool(gff_file).intersect(bed, u=True)
    # Save fasta file
    fasta_name = os.path.join(fasta_folder, str(sequence) + "_" + os.path.basename(output) + ".fa")
    fasta.save_seqs(fasta_name)
    # Save gff file
    gff_name = os.path.join(gff_folder, str(sequence) + "_" + os.path.basename(output) + ".gff")
    gff.saveas(gff_name)

    return


def extract(shortlist, start, finish, gpa, core, output, fasta_list, gff_list):
    # Extract accessory genes between core segments from ordered and annotated gene_presence_absence files
    # Check first gene is not Sequence_break which is non-unique
    if start != "Sequence_break":
        if finish == "Sequence_break":
            for sequence in shortlist.keys():
                # Read in segment
                segment = pd.read_csv(gpa + "/" + sequence + "_gene_presence_absence.csv").fillna("")
                segment = segment.rename(columns={"Non-unique Gene name": "Non-unique_name", "Annotation": "Annotation",
                                                  "No. isolates": "No_isolates", "No. sequences": "No_sequences"})
                index = segment.index[segment["Gene"] == start][0]
                # Need to add one as pandas slice excludes end index
                n_acc = shortlist[sequence]["Core_region_accessory_count"]
                index_fow = index + n_acc + 1
                index_rev = index - n_acc
                # Check to make sure index_rev is not negative - may occur if merged genes
                if index_rev < 0:
                    index_rev = 0
                # Check we will slice in right direction - sometimes flipped
                # Find contig in forward direction and check number of accessory genes in that direction
                # and how many lie on same contig as core gene
                seg_forward = segment.iloc[(index + 1):index_fow]
                chrom = segment.iloc[index, 9]
                for_acc = len(seg_forward.iloc[:, 6][(seg_forward.iloc[:, 6] < core) &
                                                     (seg_forward.iloc[:, 9] == chrom)])
                # Check reverse
                seg_rev = segment.iloc[index_rev:index]
                rev_acc = len(seg_rev.iloc[:, 6][(seg_rev.iloc[:, 6] < core) &
                                                 (seg_rev.iloc[:, 9] == chrom)])
                # Pick the direction with the most valid accessory genes
                if for_acc >= rev_acc:
                    acc_segment = seg_forward[(seg_forward["No_isolates"] < core) &
                                              (seg_forward["chromosome"] == chrom)]
                    if len(acc_segment) > 0:
                        # Set bedtools bed interval from end of first core gene to end of last accessory gene
                        bedtools_start = segment.iloc[
                            segment.index[segment["Gene"] == acc_segment.iloc[0, 0]][0] - 1, 11]
                        bedtools_end = acc_segment.iloc[-1, 11]
                        # Get fasta
                        # As segment contains a sequence break, extract regardless of presence of MGE genes
                        getfasta(bedtools_start, bedtools_end, chrom, sequence, output, fasta_list, gff_list)
                else:
                    acc_segment = seg_rev[(seg_rev["No_isolates"] < core) &
                                          (seg_rev["chromosome"] == chrom)]
                    if len(acc_segment) > 0:
                        # Set bedtools bed interval from start of first accessory gene to start of core gene
                        bedtools_start = acc_segment.iloc[0, 10] - 1
                        bedtools_end = segment.iloc[
                                           segment.index[segment["Gene"] == acc_segment.iloc[-1, 0]][0] + 1, 10] - 1
                        # Get fasta
                        # As segment contains a sequence break, extract regardless of presence of MGE genes
                        getfasta(bedtools_start, bedtools_end, chrom, sequence, output, fasta_list, gff_list)

                # Write segments to csv
                gff = segment.columns[8]
                acc_segment.to_csv(output + "/" + gff + ".filtered.temp", index=False)

        # If both genes listed
        else:
            for sequence in shortlist.keys():
                segment = pd.read_csv(gpa + "/" + sequence + "_gene_presence_absence.csv").fillna("")
                segment = segment.rename(columns={"Non-unique Gene name": "Non-unique_name", "Annotation": "Annotation",
                                                  "No. isolates": "No_isolates", "No. sequences": "No_sequences"})
                index_start = segment.index[segment["Gene"] == start][0]
                index_finish = segment.index[segment["Gene"] == finish][0]
                chrom = segment.iloc[index_start, 9]
                # Check start and end genes are on same chromosome
                if segment.iloc[index_start, 9] != segment.iloc[index_finish, 9]:
                    print("Error: " + start + " and " + end + " are not on the same chromosome for " + sequence)
                    exit(1)
                # Check which direction to slice in
                if index_start < index_finish:
                    acc_segment = segment.iloc[(index_start + 1):index_finish, :]
                else:
                    acc_segment = segment.iloc[(index_finish + 1):index_start, :]

                # Remove core gene read-through
                acc_segment = acc_segment[acc_segment["No_isolates"] < core]

                # Set bed interval
                if len(acc_segment) > 0:
                    # Set interval from gene before first acc gene and gene after last acc gene
                    # This accounts from removed core genes
                    bedtools_start = segment.iloc[segment.index[segment["Gene"] == acc_segment.iloc[0, 0]][0] - 1, 11]
                    bedtools_end = segment.iloc[
                                       segment.index[segment["Gene"] == acc_segment.iloc[-1, 0]][0] + 1, 10] - 1
                    # Get fasta
                    # Check if there's a recombinase gene in the segment
                    if len([x for x in acc_segment["recombinase"].to_list() if x != '']) >= 1:
                        getfasta(bedtools_start, bedtools_end, chrom, sequence, output, fasta_list, gff_list)
                # Write segments to csv
                gff = segment.columns[8]
                acc_segment.to_csv(output + "/" + gff + ".filtered.temp", index=False)

    # In the case of first gene being "Sequence_break", index from second gene
    else:
        for sequence in shortlist.keys():
            # Read in gene_presence_absence file
            segment = pd.read_csv(gpa + "/" + sequence + "_gene_presence_absence.csv").fillna("")
            segment = segment.rename(columns={"Non-unique Gene name": "Non-unique_name", "Annotation": "Annotation",
                                              "No. isolates": "No_isolates", "No. sequences": "No_sequences"})
            index = segment.index[segment["Gene"] == finish][0]
            # Need to add one as pandas slice excludes end index
            n_acc = shortlist[sequence]["Core_region_accessory_count"]
            index_fow = index + n_acc + 1
            index_rev = index - n_acc
            # Check to make sure index_rev is not negative - may occur if merged genes
            if index_rev < 0:
                index_rev = 0

            # Check we will slice in right direction - sometimes flipped
            # Find contig in forward direction and check number of accessory genes in that direction
            # and how many lie on same contig as core gene
            seg_forward = segment.iloc[(index + 1):index_fow]
            chrom = segment.iloc[index, 9]
            for_acc = len(seg_forward.iloc[:, 6][(seg_forward.iloc[:, 6] < core) &
                                                 (seg_forward.iloc[:, 9] == chrom)])
            # Check reverse
            seg_rev = segment.iloc[index_rev:index]
            rev_acc = len(seg_rev.iloc[:, 6][(seg_rev.iloc[:, 6] < core) &
                                             (seg_rev.iloc[:, 9] == chrom)])
            # Pick the direction with the most valid accessory genes
            if rev_acc >= for_acc:
                acc_segment = seg_rev[(seg_rev["No_isolates"] < core) &
                                      (seg_rev["chromosome"] == chrom)]
                if len(acc_segment) > 0:
                    # Set bedtools bed interval from start of first accessory gene to start of core gene
                    bedtools_start = acc_segment.iloc[0, 10] - 1
                    bedtools_end = segment.iloc[
                                       segment.index[segment["Gene"] == acc_segment.iloc[-1, 0]][0] + 1, 10] - 1
                    # Get fasta
                    # As segment contains sequence break, extract regardless of presence of MGE genes
                    getfasta(bedtools_start, bedtools_end, chrom, sequence, output, fasta_list, gff_list)
            else:
                acc_segment = seg_forward[(seg_forward["No_isolates"] < core) &
                                          (seg_forward["chromosome"] == chrom)]
                if len(acc_segment) > 0:
                    # Set bedtools bed interval from end of core gene to end of last accessory gene
                    bedtools_start = segment.iloc[segment.index[segment["Gene"] == acc_segment.iloc[0, 0]][0] - 1, 11]
                    bedtools_end = acc_segment.iloc[-1, 11]
                    # Get fasta
                    # As segment contains sequence break, extract regardless of presence of MGE genes
                    getfasta(bedtools_start, bedtools_end, chrom, sequence, output, fasta_list, gff_list)

            # Write segments to csv
            gff = segment.columns[8]
            acc_segment.to_csv(output + "/" + gff + ".filtered.temp", index=False)


def main():
    description = "Pulls accessory segment between specified core-core pair"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-f",
        "--first",
        dest="f_gene",
        help="First gene in core-core pair",
        required=True,
        type=str,
    )

    io_opts.add_argument(
        "-s",
        "--second",
        dest="s_gene",
        help="Second gene in core-core pair",
        required=True,
        type=str,
    )

    io_opts.add_argument(
        "-l",
        "--low_freq",
        dest="low_freq",
        help="low_frequency_gene_placement.tsv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-r",
        "--rec",
        dest="rec",
        help="recombinase_rules.tsv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        help="Output directory",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-fa",
        "--fasta",
        dest="fasta_dir",
        help="Directory to fasta files",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-gf",
        "--gff",
        dest="gff_dir",
        help="Directory to gff files",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-i",
        "--seq_id",
        dest="seq_id",
        help="sequence_id.txt",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-g",
        "--gpa",
        dest="gpa",
        help="Path to ordered_gene_presence_absence folder",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-c",
        "--core",
        dest="core",
        help="Core threshold expressed as fraction up to 1",
        required=False,
        type=float,
        default=0.99
    )

    io_opts.add_argument(
        "-m",
        "--min",
        dest="min",
        help="Minimum accessory content",
        required=False,
        type=int,
        default=1
    )

    args = parser.parse_args()

    # Make dedicated directory for outputs
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # Make directory for core-core combination
    output = os.path.join(args.output_dir, args.f_gene + "-" + args.s_gene)
    if not os.path.exists(output):
        os.mkdir(output)

    # Make directory for fasta and gffs
    seg_fasta_folder = os.path.join(args.output_dir, "segment_fasta")
    if not os.path.exists(seg_fasta_folder):
        os.mkdir(seg_fasta_folder)
    # Make folder for gff files
    seg_gff_folder = os.path.join(args.output_dir, "segment_gff")
    if not os.path.exists(seg_gff_folder):
        os.mkdir(seg_gff_folder)

    # read in low_frequency_gene_placement.tsv
    low_freq = pd.read_csv(args.low_freq, sep="\t")
    low_freq = low_freq[(low_freq["Core_gene_1"] == args.f_gene) &
                        (low_freq["Core_gene_2"] == args.s_gene) &
                        (low_freq["Core_region_accessory_count"] >= args.min)]

    # Create an ordered sequence list
    # This will be merged later to create a master list
    seq_id = pd.read_csv(args.seq_id, header=None, names=["Gff"])
    # Left merge low_freq with seq_id to maintain order of Gffs
    joined = seq_id.merge(low_freq, how="left", on="Gff")

    # Change to integer
    joined[["Core_region_accessory_count", "Core_region_size"]] = \
        joined[["Core_region_accessory_count", "Core_region_size"]].fillna(0).astype(int)

    # Output accessory gene count
    acc_joined = joined[["Gff", "Core_region_accessory_count"]].astype(str).replace("0", "")
    acc_joined = acc_joined.rename({"Core_region_accessory_count": (args.f_gene + "-" + args.s_gene)}, axis=1)
    acc_joined.to_csv(output + "/" + args.f_gene + "-" + args.s_gene + "_joined_count.tsv", sep="\t", index=False)
    dist_joined = joined[["Gff", "Core_region_size"]].astype(str).replace("0", "")
    dist_joined = dist_joined.rename({"Core_region_size": (args.f_gene + "-" + args.s_gene)}, axis=1)
    dist_joined.to_csv(output + "/" + args.f_gene + "-" + args.s_gene + "_joined_distance.tsv", sep="\t", index=False)

    # Calculate number of core genes and core gene cutoff
    core = np.floor(len(seq_id.index) * args.core)

    # Sequences to examine
    shortlist = low_freq[["Gff", "Core_region_accessory_count"]].set_index("Gff").to_dict("index")

    # List of fasta files - searches for .fa and .fasta extensions
    fasta_list = [os.path.join(args.fasta_dir, x) for x in os.listdir(args.fasta_dir) if
                  x.endswith('.fa') or x.endswith('.fasta')]

    # List of gff files - searches for .fa and .fasta extensions
    gff_list = [os.path.join(args.gff_dir, x) for x in os.listdir(args.gff_dir) if x.endswith('.gff')]

    # Extract accessory genes between core segments from ordered and annotated gene_presence_absence files
    extract(shortlist, args.f_gene, args.s_gene, args.gpa, core, output, fasta_list, gff_list)


if __name__ == "__main__":
    main()
