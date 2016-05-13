#!/usr/local/bin/python2.7

#####################
##
## Imports
##

import sys
import re
import argparse
import numpy
import copy
from chopchop.data import connect_db, geneToCoord_db, geneToCoord_BED, regionsToFasta, print_fasta, read_json_file, convert_to_relative_coordinates, print_genbank, print_bed


def formatDNA(dnaSeq, start, coords):
    
    coords = numpy.array(coords)
    
    # Create array of coordinates for exons, zero indexed
    zeroIndexCoords = coords - start
    
    formatSeq = ""
    
    # Check whether each base lies within the exon coordinates. If yes, convert to uppercase
    for base in range(0, len(dnaSeq)):
        for exon in zeroIndexCoords:
        
            if base in range(exon[0], exon[1]):            
                formatSeq+=dnaSeq[base].upper()
                break
                
        else:
            formatSeq+=dnaSeq[base]
                
    return formatSeq
    
    


## Recursive merge
def merge(el):
    if len(el) == 1:
        return el

    if el[0][2] > el[1][1]:
        newEl = [(el[0][0], el[0][1], el[1][2])]
        newEl.extend(el[2:len(el)])

        return merge(newEl)
    else:
        newEl = [el[0]]
        newEl.extend(merge(el[1:len(el)]))
        return newEl



def add_flanks(regions, flanks):
    regions = [(x[0], int(x[1])-flanks, int(x[2])+flanks) for x in regions]

    # Merge regions that overlap
    return merge(regions)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", help="Target genes or regions", metavar="TARGET_REGIONS")
    parser.add_argument("twoBitFile", help="File containing the 2bit compressed genome.", metavar="TWOBIT_FILE")
    parser.add_argument("outputFile", help="File to print the output.", metavar="OUTPUT_FILE")
    parser.add_argument("-O", "--output", help="Type of output file: FASTA/GENBANK/BED/GTF.", default="FASTA")   
    parser.add_argument("-i", "--gene_index", help="If there are multiple entries with the same name use number INDEX.", default=1, type=int, metavar="INDEX")   
    parser.add_argument("-I", "--include_introns", help="Include intron regions", default=False, action="store_true")   
    parser.add_argument("-f", "--flanks", help="Include flanking regions", type=int, default=0)   
    parser.add_argument("-F", "--fasta", default=False, action="store_true", help="Use FASTA file as input rather than gene or genomic region.")
    parser.add_argument("-G", "--genome", default="danRer7", metavar="GENOME", help="The genome to search. Used for gene index.") 
    parser.add_argument("-g", "--gene_file", help="Use a BED file to retrieve gene information.", metavar="BED_FILE", dest="bed_file")   
    parser.add_argument("-D", "--database", help="Connect to a chopchop database to retrieve gene: user_name:passwd@host/database", metavar="DATABASE", dest="database")   
    parser.add_argument("-d", "--description", help="Adds a description to the output file", metavar="DESCRIPTION", default="")   
    parser.add_argument("-P", "--primers", help="The file specifying the coordinates of the primers.", metavar="PRIMER_FILE", dest="primers")   
    parser.add_argument("-T", "--targetAnnotation", help="The file specifying the coordinates of the targets.", metavar="TARGET_FILE", dest="targetAnnotation")   
    parser.add_argument("-V", "--targetVis", help="The file specifying the coordinates of the targets.", metavar="TARGET_FILE", dest="targetVis")  
    parser.add_argument("-m", "--mode", help="Mode.", default="CRISPR", dest="mode")

    args = parser.parse_args()

    # Connect to database if requested
    if args.database:
        cdb = connect_db(args.database)
        db = cdb.cursor()
        use_db = True
    else:
        use_db = False


    if args.fasta:
        sequence = ""
        fastaFile = open(args.targets, 'r')
        for line in fastaFile:
            if line[0] == ">":
                seq_name = line[1:].strip()
            else:
                sequence += line.rstrip()        
        fastaFile.close()

        targets = [["sequence", 0, len(sequence)]]
        
        regions = targets
        targetRegions = None
    else:
        
    ## Check if coordinate
        pattern = re.compile("((\w+):)?([\.\,\d]+)\-([\.\,\d]+)")
        isCoordinate = pattern.match(args.targets)

        if isCoordinate: 
            coords = [[isCoordinate.group(2), int(isCoordinate.group(3)), int(isCoordinate.group(4))]]
            regions = coords
            targetRegions = copy.deepcopy(regions)
#        sys.stderr.write("\n%s\n\n" % ",".join(coords))
#        pass
        else:
            if use_db:
                coords = geneToCoord_db(args.targets, args.genome, db, args.gene_index)
            else:
                if args.bed_file:
                    coords = geneToCoord_BED(args.targets, args.bed_file, args.gene_index)
                else:
                    sys.stderr.write("Need to provide a source, either file (-g) or database (-D).\n")
                    sys.exit(1)


        # Get regions for annotation
            startSplit = coords[1].split(",")
            endSplit = coords[2].split(",")
            del startSplit[-1]
            del endSplit[-1]                
            startSplit = map(int, startSplit)
            endSplit = map(int, endSplit)
            targetRegions = zip([coords[0]] * len(startSplit),startSplit,endSplit)
            regions = targetRegions

        # Return whole locus if include introns
            if args.include_introns:
                regions = [(coords[0], coords[7], coords[8])]

        ## ADD flanks    
        regions = add_flanks(regions, args.flanks)


        # Get FASTA from 2bit index
        #    sys.stderr.write("BEFORE\n")
        sequence = regionsToFasta(regions, args.twoBitFile, targetRegions)
        #    sys.stderr.write("AFTER\n")
        #    sys.exit()


    if args.output == "FASTA":
        print_fasta(args.targets, sequence, args.outputFile, args.description)
    elif args.output == "GENBANK":

        targets = read_json_file(args.targetAnnotation)
        targets = convert_to_relative_coordinates(regions, targets)

        if targetRegions != None:
            exons = [list(elem) for elem in targetRegions]
            exons = convert_to_relative_coordinates(regions, exons, endCoord=2)
        else:
            exons = None

        print_genbank(args.targets, sequence, exons, None, targets, args.outputFile, args.description)
    elif args.output == "BED":

        targets = read_json_file(args.targetAnnotation)
        vis = read_json_file(args.targetVis)
        print_bed(args.mode, args.targets, vis, targets, args.outputFile, args.description)


if __name__ == '__main__':
    main()
