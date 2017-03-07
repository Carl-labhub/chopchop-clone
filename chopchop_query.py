#!/usr/bin/env python2.7

import os
import csv
import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--genePred_file",
                       help="Full path to genePred file from which gene names will be extracted.")
    group.add_argument("--gene_names",
                       help="You can supply comma separated list of gene names eg. rho,mtm2.")
    parser.add_argument("-o", "--outputDir",
                        help="Output directory. It is important to not run this script without creating beforehand a folder for results. Cleanup after CHOPCHOP scripts will remove everything with the exception of .tsv files.")
    args = parser.parse_known_args()

    if args[0].gene_names is not None:
        genes = args[0].gene_names.split(",")
    else:
        genes = []
        tableR = open(args[0].genePred_file, 'rb')
        tablereader = csv.DictReader(
            tableR, delimiter='\t', quoting=csv.QUOTE_NONE)

        # Look in genome table for gene of question
        for row in tablereader:
            genes.append(row['name'])

        tableR.close()

    if args[0].outputDir == "./":
        exit("Supply a folder for results.")

    relative_cc = os.path.join(os.path.dirname(__file__), "chopchop.py")
    command = args[1]
    command.insert(0, relative_cc)
    command.append("-o")
    command.append(args[0].outputDir)

    for gene in genes:

        gene_file = os.path.join(args[0].outputDir, gene + ".tsv")
        gene_handle = open(gene_file, "w")

        temp_comd = list(command)
        temp_comd.append(gene)

        b = subprocess.Popen(temp_comd,
                             stdout=gene_handle, stderr=subprocess.STDOUT)
        b.wait()

        print gene

    # cleanup after CHOPCHOP temp files
    for the_file in os.listdir(args[0].outputDir):
        file_path = os.path.join(args[0].outputDir, the_file)
        try:
            if os.path.isfile(file_path) and os.path.splitext(the_file)[1] != ".tsv":
                os.unlink(file_path)
        except Exception as e:
            print(e)


if __name__ == '__main__':
    main()
