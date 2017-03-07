#!/usr/bin/env python2.7

#####################
##
## Imports
##
import re
import os
import sys
import csv
import copy
import math
import json
import string
import argparse
import pickle
import pandas
import numpy
import featurization as feat
import scipy.stats as ss

from collections import defaultdict
from Bio import SeqIO, SeqFeature, SeqRecord
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import FeatureLocation
from operator import itemgetter, attrgetter
from subprocess import Popen, PIPE


#####################
##
## Global variables
##
ISOFORMS = False

# PATHs
PRIMER3 = "./primer3_core"
BOWTIE = "bowtie/bowtie"
TWOBITTOFA = "./twoBitToFa"
TWOBIT_INDEX_DIR = "/home/ai/Projects/data/genomes_chopchop"
BOWTIE_INDEX_DIR = "/home/ai/Projects/data/genomes_chopchop"
ISOFORMS_INDEX_DIR = "/home/ai/Projects/c2c2"
GENE_TABLE_INDEX_DIR = "/home/ai/Projects/c2c2"

# Program mode
CRISPR = 1
TALENS = 2
CPF1 = 3
NICKASE = 4

# Maximum genomic region that can be searched
TARGET_MAX = 20000

# Defaults
CRISPR_DEFAULT = {"GUIDE_SIZE" : 20,
                  "PAM": "NGG",
                  "MAX_OFFTARGETS" : 50,
                  "MAX_MISMATCHES" : 3,
                  "SCORE_GC" : True,
                  "SCORE_FOLDING" : True}

TALEN_DEFAULT = {"GUIDE_SIZE" : 18,
                 "PAM": "",
                 "MAX_OFFTARGETS" : 100,
                 "MAX_MISMATCHES" : 2,
                 "SCORE_GC" : False,
                 "SCORE_FOLDING" : False}

CPF1_DEFAULT =  {"GUIDE_SIZE" : 24,
                 "PAM": "TTTN",
                 "MAX_OFFTARGETS" : 50,
                 "MAX_MISMATCHES" : 3,
                 "SCORE_GC" : True,
                 "SCORE_FOLDING" : True}

NICKASE_DEFAULT =  {"GUIDE_SIZE" : 20,
                    "PAM": "NGG",
                    "MAX_OFFTARGETS" : 50,
                    "MAX_MISMATCHES" : 3,
                    "SCORE_GC" : True,
                    "SCORE_FOLDING" : True}

TALEN_OFF_TARGET_MIN = 28
TALEN_OFF_TARGET_MAX = 42
PRIMER_OFF_TARGET_MIN = 1
PRIMER_OFF_TARGET_MAX = 1000

# Max members of a TALENs cluster (15)
MAX_IN_CLUSTER = 15

# SCORES
SCORE = {"INPAIR_OFFTARGET_0" : 5000,
         "INPAIR_OFFTARGET_1" : 3000,
         "INPAIR_OFFTARGET_2" : 2000,
         "INPAIR_OFFTARGET_3" : 1000,
         "OFFTARGET_PAIR_SAME_STRAND" : 10000,
         "OFFTARGET_PAIR_DIFF_STRAND" : 5000,
         "MAX_OFFTARGETS" : 20000, ## FIX: SPECIFIC FOR TALEN AND CRISPR
         "COEFFICIENTS" : 100,
         "CRISPR_BAD_GC" : 500,
         "FOLDING" : 1}

SINGLE_OFFTARGET_SCORE = [1000, 800, 600, 400]
GC_LOW = 40
GC_HIGH = 70

G_20 = {"Intercept": -30,
        "G1": 60}

XU_2015 = {'C18':-0.113781378,
           'G17':0.080289971,
           'A16':0.025840846,'G16':0.072680697,
           'G15':0.100642827,
           'G14':0.082839514,
           'T14':-0.070933894,
           'A12':0.02156311,
           'A11':0.129118902,
           'A10':0.030483786,'T10':-0.169986128,
           'A9':0.093646913,
           'G7':-0.214271553,'T7':0.073750154,
           'A6':0.202820147,
           'A5':0.129158071,
           'G4':0.107523301,'T4':-0.349240474,
           'C3':0.23502822,'T3':-0.145493093,
           'G2':0.238517854,'T2':-0.300975354,
           'C1':-0.125927965,'G1':0.353047311,'T1':-0.221752041,
           'PAMT1':-0.155910373,
           '1C':0.179639101,
           '4T':-0.116646129}

DOENCH_2014 = {"Intercept": 0.5976361543,
               "G23": -0.2753771278,"TG22": -0.625778696,
               "A22": -0.3238874564,"C22": 0.1721288713,
               "C21": -0.1006662089,
               "C20": -0.20180294,"G20": 0.2459566331,"CG19": 0.3000433167,
               "C19": 0.0983768352,"A19": 0.0364400412,"AA18":-0.8348362447,"AT18": 0.7606277721,
               "C18":-0.7411812913,"G18":-0.3932643973,"GG17":-0.4908167494,
               "A13":-0.4660990147,"GG12":-1.5169074394,"AT12": 0.7092612002,"CT12": 0.4962986088,"TT12":-0.5868738941,
               "GG11":-0.3345637351,
               "AG10": 0.7638499303,"CG10":-0.5370251697,
               "A10": 0.0853769455,"C10":-0.0138139718,
               "A9": 0.2726205124,"C9": -0.119022648,"T9":-0.2859442224,
               "A8": 0.0974545916,"G8":-0.1755461698,"GT7":-0.7981461328,
               "C7":-0.3457954508,"G7":-0.6780964263,
               "A6": 0.2250890296,"C6":-0.5077940514,"GG5":-0.6668087295,"CT5": 0.3531832525,
               "G5":-0.4173735974,"T5":-0.0543069593,"CC4": 0.7480720923,"GT4":-0.3672667722,
               "G4": 0.379899366,"T4":-0.0907126437,"CA3": 0.5682091316,"GC3": 0.3290720742,"AG3":-0.8364567552,"GG3":-0.7822075841,
               "C3": 0.0578233185,"T3":-0.5305672958,"CT2":-1.0296929571,
               "T2":-0.8770074285,"GC1": 0.8561978226,"TC1":-0.4632076791,
               "C1":-0.8762358461,"G1": 0.2789162593,"T1":-0.4031022177,"AA0":-0.5794923887,"GA0": 0.6490755373,
               "PAMC1": 0.287935617,"PAMA1":-0.0773007042,"PAMT1":-0.2216372166, "PAMAG1":-0.0773007042,"PAMCG1":  0.287935617,"PAMTG1":-0.2216372166,
               "1G":-0.6890166818,"1T": 0.1178775773,
               "2C":-0.1604453039,"2GG":-0.6977400239,
               "3G": 0.3863425849,
               "gc_low":-0.2026258943,
               "gc_high": -0.166587752}

MORENO_MATEOS_2015 = {"Intercept": 0.1839309436,
                      "G26":-0.0296937089,
                      "CG23":0.0246817853,"GT23":0.0229499956,
                      "G23":-0.0054488693,"A23":-0.0421535206,
                      "C18":0.0024492239,"G18":0.1146006812,"GG17":-0.0015779899,"CG17":0.0541714023,
                      "G17":0.0677391822,"GA16":0.0637170933,"GG16":0.0268021579,"AA16":-0.0169054146,
                      "A16":-0.0182872921,"G16":0.0209290394,"TG15":0.0536784362,
                      "A15":0.0116332345,"G15":0.0275911379,"GG14":0.0418830086,
                      "A14":0.0176289243,"T14":0.0354451707,"C14":0.069495944,"G14":0.0613609047,"GG13":0.0743558476,"TT13":-0.0861877104,
                      "G13":0.0251167144,"A13":-0.0184872292,
                      "A12":-0.0105952955,"C12":-0.0004777273,"G12":0.0511297167,"GT11":0.0533728222,
                      "G11":0.0379709424,"C11":-0.0216386089,
                      "G10":0.0154937801,"TT9":0.0349288099,
                      "A9":-0.033820432,"G9":0.0164578159,"GT8":0.0459089908,"TG8":0.0023917441,"TT8":-0.094424075,
                      "A8":-0.0155764989,"G8":0.0179168437,"AA7":-0.0973770966,
                      "C7":0.0150895135,"AG6":0.0097407989,
                      "T6":-0.0687304967,"C6":0.0342629207,"CC5":0.0889196009,
                      "T5":0.0132240349,"G5":0.1011443803,"C5":0.0376316197,"A5":0.0319309088,
                      "T4":-0.0014222433,"CC3":0.0950722865,"TG3":0.1067185626,"GA3":-0.0543384557,"GT3":-0.0663880754,
                      "T3":-0.0119961724,"A3":0.0374664775,"C3":0.0529723137,"G3":0.1054883249,"AC2":0.0622193698,"TG2":0.0609521143,
                      "C2":-0.031648353,"A2":0.010506405,"GG1":0.1115594407,"CG1":-0.0734536087,
                      "G1":0.0361466487,"C1":-0.0003689729,"TC0":-0.0842648932,
                      "PAMT1":-0.0002808449,"PAMA1":0.0191268154,"PAMC1":0.0799339215,"PAMG1":0.0851510516,
                      "1G":-0.0463159143,"1C":-0.0131827326,"1T":0.0172631618,"1CA":0.0577598507,
                      "2C":-0.0307155561,"2A":0.0015897498,"2TG":0.0481368123,"2GT":0.0734253504,"2GA":-0.01227989,
                      "3G":0.0307124897,
                      "5G":-0.0141671226,"5T":-0.0176476917,"5GA":-0.0377977074,"5AG":-0.0419359085,
                      "6A":0.0485962592}

#n":umber of nucleotides counted towards scoring after PAM, end before 5' gRNA start eg. 3: 321 gRNA PAM 123
DOWNSTREAM_NUC = 15

# EXIT CODES
EXIT = {"PYTHON_ERROR" : 1,
        "BOWTIE_ERROR" : 2,
        "TWOBITTOFA_ERROR" : 3,
        "GENE_ERROR" : 4,
        "DB_ERROR" : 5,
        "PRIMER3_ERROR" : 6,
        "BOWTIE_PRIMER_ERROR" : 7,
        "ISOFORM_ERROR" : 8}


# PRIMER3 OPTIONS
PRIMER3_CONFIG = {"PRIMER_OPT_SIZE" : "22",
                  "PRIMER_MIN_SIZE" : "18",
                  "PRIMER_MAX_SIZE" : "25",
                  "PRIMER_MAX_NS_ACCEPTED" : "0",
                  "PRODUCT_SIZE_MIN" : "100",
                  "PRODUCT_SIZE_MAX" : "290"}


# WIDTH
VIS_WIDTH = {"UTR" : 10,
             "CDS" : 20}


# SELF-COMPLEMENTARITY
STEM_LEN = 4

#####################
##
## Classes
##

class Hit:
    """Creates class for each hit from bowtie."""    

    def __init__(self, line):    
        self.flagSum = int(line[1])
        self.chrom = line[2]
        self.start = int(line[3])
        self.matchSeq = line[9]
        self.mismatch = line[-1]
        self.mismatchPos = line[-2]
        self.opts = line[11:(len(line))]
        self.mismatchCorrected = False


    def calc_mismatchPos (self):        
        """ Updates the sequence parsed from the SAM output to include the mismatches """

        lastDigit = len(self.mismatchPos)-1
        guideSize = len(self.matchSeq)
        guideCurr = ""

        ## MD:Z:GUIDESIZE means that there are no mismatches
        if not(self.mismatchPos =="MD:Z:%s" % guideSize):
            guideIndex = 0
            currTotal = 0

            for c in range(5, lastDigit+1):

                # If the character is a digit, check if the next character is a digit (>9) and add number to total
                if self.mismatchPos[c].isdigit():
                    
                    if c != lastDigit and self.mismatchPos[c+1].isdigit():
                        currTotal += (int(self.mismatchPos[c])*10)
                    else:
                        currTotal += int(self.mismatchPos[c])
                        guideCurr += self.matchSeq[guideIndex:currTotal]
                        guideIndex = currTotal
                        
                # if character is a letter, add one to total
                else:
                    guideCurr += self.mismatchPos[c].lower()
                    currTotal += 1
                    guideIndex += 1
               

            self.matchSeq = guideCurr
            
  
    # Specifying how to print items in list of off-targets   
    def __str__(self):
        if not self.mismatchCorrected:
            self.calc_mismatchPos()
            self.mismatchCorrected = True
            
        return "%s:%s\t%s" % (self.chrom, self.start, self.matchSeq)

    def asOffTargetString(self, label, maxOffTargets):
        if self.mismatch == "XM:i:%s" % maxOffTargets:               
            return "%s,>%s across the genome,0-2,n/a " % (label, maxOffTargets)
        else:
            if not self.mismatchCorrected:
                self.calc_mismatchPos()
                self.mismatchCorrected = True

        return "%s,%s,%s,%s" % (label, self.chrom + ":" + str(self.start), self.mismatch[-1], self.matchSeq)


class Guide(object):
    """ This defines a class for each guide. The (off-target) hits for each guide form a separate class. The functions "addOffTarget" and
    "sort_offTargets" applies to just the Tale class """

    def __init__(self, name, flagSum, guideSize, guideSeq, scoreGC, scoreSelfComp, 
                 backbone_regions, PAM, replace5prime=None, scoringMethod=None, 
                 genome=None, gene=None, isoform=None, gene_isoforms=None):
        
        self.scoringMethod = scoringMethod
        self.genome = genome
        self.gene = gene
        self.isoform = isoform
        self.gene_isoforms = gene_isoforms
        self.offTargetsIso = {0: set(), 1: set(), 2: set(), 3: set()}
        self.conserved = False # conservation of guide across all isoforms
        self.PAM = PAM
        # From the guide's name we can get the chromosome
        self.flagSum = flagSum
        elements = name.split(":")
        self.ID = elements[0]
        self.chrom = elements[1]
        coord = elements[2]
        
        self.name = ":".join(elements[0:3])
        
        if len(elements) > 3:
            self.downstream5prim = elements[3]
            self.downstream3prim = elements[4]
            self.strand = elements[5]
        else:
            self.downstream5prim = ''
            self.downstream3prim = ''
            self.strand = None
            
        self.guideSize = guideSize
        self.targetSize = guideSize
        self.cluster = -1
        self.score = 0
        self.ALL_scores = [0, 0, 0, 0, 0, 0]

        # Off target count
        self.offTargetsMM = [0] * 4

        # The location of the last digit of the exon start in the name string
        mid = coord.find('-')        
        # The location of the first digit of the guide position in the exon
        end = (coord.find('_'))+1
                
        # The full position of the guide in the exon
        region = coord[end:]
        
        # The location of the last digit of the guide position in the exon
        location = region.find('-')

        # The start of the exon containing the guide
        self.exonStart = int(coord[0:mid])
        self.exonNum = None

        # The number of bases after the exon start
        guidePos = int(region[:location])+1

        # guide start coordinate
        self.start = self.exonStart + guidePos
        self.end = self.start + guideSize
        self.guideSeq = guideSeq

        # Record which strand the guide is on
        if flagSum == "16" or ISOFORMS: # due to reverse complementing before alignments
            self.strandedGuideSeq = guideSeq
            if self.strand is None:
                self.strand = '+'
        else:
            self.strandedGuideSeq = str(Seq(guideSeq).reverse_complement())
            if self.strand is None:
                self.strand = '-'

        # Initiate offTargets list
        self.offTargets = []
        self.offTarget_hash = {}
        self.offTargets_sorted = False

        if scoreSelfComp:
            self.calcSelfComplementarity(scoreSelfComp, backbone_regions, PAM, replace5prime)
        else:
            self.folding ="N/A"
            
        # Scoring
        self.calcGCContent(scoreGC)


    def calcSelfComplementarity(self, scoreSelfComp, backbone_regions, PAM, replace5prime = None):   
        if replace5prime:
            fwd = self.strandedGuideSeq[len(PAM):-len(replace5prime)] + replace5prime #Replace the 2 first bases with e.g. "GG"
        else:
            fwd = self.guideSeq[len(PAM):] # Do not include PAM motif in folding calculations

        rvs = str(Seq(fwd).reverse_complement())
        L = len(fwd)-STEM_LEN-1            

        self.folding = 0
        
        for i in range(0,len(fwd)-STEM_LEN):
            if gccontent(fwd[i:i+STEM_LEN]) >= 0.5:
                if fwd[i:i+STEM_LEN] in rvs[0:(L-i)] or any([fwd[i:i+STEM_LEN] in item for item in backbone_regions]):
                    #sys.stderr.write("%s\t%s\n" % (fwd, fwd[i:i+STEM_LEN]))
                    self.folding += 1
                    
        self.score += self.folding * SCORE['FOLDING']


    def calcGCContent(self, scoreGC):
        """ Calculate the GC content of the guide """
        Gcount = self.guideSeq.count('G')
        Ccount = self.guideSeq.count('C')
        self.GCcontent = (100*(float(Gcount+Ccount)/int(self.guideSize)))

        if scoreGC:
            if self.GCcontent > GC_HIGH or self.GCcontent < GC_LOW:
                self.score += SCORE['CRISPR_BAD_GC']


    def addOffTarget(self, hit, maxOffTargets, countMMPos):        
        """ Add off target hits (and not original hit) to list for each guide RNA """    

        hit_id = "%s:%s" % (hit.chrom, hit.start)
        nmiss = 0
        mm_pattern = re.compile('NM:i:(\d+)')
       
        # If the hit is identical to the guide coord it is the original correct hit
        if self.chrom == hit.chrom and self.start == hit.start: # never true for isoforms
            # This is the original/main hit
            self.correct_hit = hit
            return
        
        if ISOFORMS and self.isoform == hit.chrom and self.strandedGuideSeq == hit.matchSeq:
            # This is the original/main hit
            self.correct_hit = hit
            return

        # Do not count off targets twice, e.g. for TALENs valid on both strands.
        if self.offTarget_hash.has_key(hit_id):  
            return

        # Reverse count+allowed arrays if on the reverse strand
        if hit.flagSum == 0 and not ISOFORMS:
            countMMPos = countMMPos[::-1]

        self.offTarget_hash[hit_id] = hit         
        MMs = get_mismatch_pos(hit.mismatchPos[5:])

        for mm in MMs:
            if not countMMPos[mm]:
                del(self.offTarget_hash[hit_id])
                return

            elif not countMMPos[mm]:
                nmiss += 1

        # Calculate score
        for opt in hit.opts:
            m = mm_pattern.match(opt)
            if m:
                mm = int(m.group(1)) - nmiss
                
                # ugly repeat to save time from iterating all isoforms 
                if ISOFORMS: 
                    if hit.chrom in self.gene_isoforms: # and hit.chrom not in self.offTargetsIso[mm]:
                        self.offTargetsIso[mm].add(hit.chrom)
                        # don't count/score isoform mismatches but display which isoforms have them
                    else: 
                        self.offTargetsMM[mm] += 1
                        self.score += SINGLE_OFFTARGET_SCORE[mm]
                else: 
                    self.offTargetsMM[mm] += 1
                    self.score += SINGLE_OFFTARGET_SCORE[mm]
                
            if opt == "XM:i:" + str(maxOffTargets):  
                self.score += SCORE['MAX_OFFTARGETS']                
                self.offTargetsMM[0] += maxOffTargets
                self.offTargetsMM[1] += maxOffTargets
                self.offTargetsMM[2] += maxOffTargets
                self.offTargetsMM[3] += maxOffTargets

        self.offTargets_sorted = False


    def numOffTargets(self):
        """ Returns the number of off-target hits for each guide """
        self.sort_offTargets()
        return len(self.offTargets)

        
    def sort_offTargets(self):
        """ Sort off-target hits according to chromosome and genomic coordinate """

        if self.offTargets_sorted:
            return        

        self.offTargets = self.offTarget_hash.values()
        self.offTargets = sorted(self.offTargets, key=attrgetter('chrom', 'start'))  
        self.offTargets_sorted = True


    def __str__(self):
        self.sort_offTargets()
        if ISOFORMS:
            return "%s\t%s:%s\t%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.strandedGuideSeq, 
                                                                                            self.chrom, self.start, self.exonNum, 
                                                                                            self.gene, self.isoform,
                                                                                            self.GCcontent, self.folding, 
                                                                                            self.offTargetsMM[0], self.offTargetsMM[1], 
                                                                                            self.offTargetsMM[2], self.offTargetsMM[3], 
                                                                                            self.conserved, (",").join(set(self.offTargetsIso[0])), 
                                                                                            (",").join(set(self.offTargetsIso[1])), 
                                                                                            (",").join(set(self.offTargetsIso[2])), 
                                                                                            (",").join(set(self.offTargetsIso[3])))
        return "%s\t%s:%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s" % (self.strandedGuideSeq, self.chrom, self.start, self.exonNum, 
                                                                self.strand, self.GCcontent, self.folding, 
                                                                self.offTargetsMM[0], self.offTargetsMM[1], 
                                                                self.offTargetsMM[2], self.offTargetsMM[3])
                  

    def asOffTargetString(self, label, maxOffTargets):
        self.sort_offTargets()
        offTargets = map(lambda x: x.asOffTargetString(label, maxOffTargets), self.offTargets)

        return ";".join(offTargets)


class Cas9(Guide):
    def __init__(self, *args, **kwargs):
        super(Cas9, self).__init__(*args, **kwargs)
        self.CoefficientsScore = {"XU_2015": 0,
                                  "DOENCH_2014": 0,
                                  "DOENCH_2016": 0,
                                  "MORENO_MATEOS_2015": 0,
                                  "CHARI_2015": 0,
                                  "G_20": 0}

        if self.scoringMethod in ["DOENCH_2016", "ALL"]:
            self.CoefficientsScore["DOENCH_2016"] = scoreDoench_2016(self.downstream5prim + self.strandedGuideSeq[:-len(self.PAM)], self.strandedGuideSeq[-len(self.PAM):], self.downstream3prim)
            if self.scoringMethod == "DOENCH_2016":
                self.score = self.score - self.CoefficientsScore["DOENCH_2016"] * SCORE['COEFFICIENTS']

        if self.scoringMethod not in ["CHARI_2015", "DOENCH_2016", "ALL"]:
            self.CoefficientsScore[self.scoringMethod] = scoregRNA(self.downstream5prim + self.strandedGuideSeq[:-len(self.PAM)], self.strandedGuideSeq[-len(self.PAM):], self.downstream3prim, globals()[self.scoringMethod])
            self.score = self.score - self.CoefficientsScore[self.scoringMethod] * SCORE['COEFFICIENTS']

        if self.scoringMethod == "ALL":
            for met in ["XU_2015", "DOENCH_2014", "MORENO_MATEOS_2015", "G_20"]:
                self.CoefficientsScore[met] = scoregRNA(self.downstream5prim + self.strandedGuideSeq[:-len(self.PAM)], self.strandedGuideSeq[-len(self.PAM):], self.downstream3prim, globals()[met])
    
    def __str__(self):
        self.sort_offTargets()
        if self.scoringMethod == "ALL":
            return "%s\t%s:%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (self.strandedGuideSeq, 
                                                                                                        self.chrom, 
                                                                                                        self.start, 
                                                                                                        self.exonNum, 
                                                                                                        self.strand, 
                                                                                                        self.GCcontent, 
                                                                                                        self.folding, 
                                                                                                        self.offTargetsMM[0], 
                                                                                                        self.offTargetsMM[1], 
                                                                                                        self.offTargetsMM[2], 
                                                                                                        self.offTargetsMM[3], 
                                                                                                        self.CoefficientsScore["XU_2015"],
                                                                                                        self.CoefficientsScore["DOENCH_2014"],
                                                                                                        self.CoefficientsScore["DOENCH_2016"],
                                                                                                        self.CoefficientsScore["MORENO_MATEOS_2015"],
                                                                                                        self.CoefficientsScore["CHARI_2015"],
                                                                                                        self.CoefficientsScore["G_20"])
        else:
            return "%s\t%s:%s\t%s\t%s\t%.0f\t%s\t%s\t%s\t%s\t%s\t%.2f" % (self.strandedGuideSeq, 
                                                                          self.chrom, 
                                                                          self.start, 
                                                                          self.exonNum, 
                                                                          self.strand, 
                                                                          self.GCcontent, 
                                                                          self.folding, 
                                                                          self.offTargetsMM[0], 
                                                                          self.offTargetsMM[1], 
                                                                          self.offTargetsMM[2], 
                                                                          self.offTargetsMM[3], 
                                                                          self.CoefficientsScore[self.scoringMethod])
    
    def calcSelfComplementarity(self, scoreSelfComp, backbone_regions, PAM, replace5prime = None):
        if replace5prime:
            fwd = replace5prime + self.strandedGuideSeq[len(replace5prime):-len(PAM)] # Replace the 2 first bases with e.g. "GG"
        else:
            fwd = self.guideSeq[0:-len(PAM)] # Do not include PAM motif in folding calculations

        rvs = str(Seq(fwd).reverse_complement())
        L = len(fwd)-STEM_LEN-1            

        self.folding = 0
        
        for i in range(0,len(fwd)-STEM_LEN):
            if gccontent(fwd[i:i+STEM_LEN]) >= 0.5:
                if fwd[i:i+STEM_LEN] in rvs[0:(L-i)] or any([fwd[i:i+STEM_LEN] in item for item in backbone_regions]):
                    #sys.stderr.write("%s\t%s\n" % (fwd, fwd[i:i+STEM_LEN]))
                    self.folding += 1
                    
        self.score += self.folding * SCORE['FOLDING']
        
        
class Pair:
    """ Pair class for 2 TALEs that are the correct distance apart """
    def __init__(self, tale1, tale2, spacerSeq, spacerSize, offTargetPairs, enzymeCo, maxOffTargets, g_RVD, minResSiteLen):
        self.tale1 = tale1
        self.tale2 = tale2
        self.chrom = tale1.chrom 
        self.strand = tale1.strand
        self.ID = ""
        self.tale1.rvd = ""
        self.tale2.rvd = ""
        self.restrictionSites = ""

        # Start of region covered by tale pair
        self.start = tale1.start

        # End of region covered by tale pair
        self.end = tale2.end # + tale2.guideSize
        self.spacerSeq = spacerSeq
        self.targetSize = spacerSize
        self.spacerSize = spacerSize
        self.offTargetPairs = offTargetPairs
        self.diffStrandOffTarget = 0
        self.sameStrandOffTarget = 0        


        # Start cluster as -1, but will increment from 1
        self.cluster = -1
        self.spacerStart = tale1.start + tale1.guideSize
        self.spacerEnd = tale2.start - 1

        self.enzymeCo = enzymeCo
        self.strandedGuideSeq = str(self.tale1.guideSeq) + "\n" + self.spacerSeq + "\n" + str(self.tale2.guideSeq)

        # Calculate RVD for TALEs; FIX: use mapping
        for base in tale1.guideSeq:
            if base == "A":
                tale1.rvd += "NI "
            elif base == "T":
                tale1.rvd += "NG "
            elif base == "C":
                tale1.rvd += "HD "
            elif base == "G":
                tale1.rvd += g_RVD
                
        for base in Seq(tale2.guideSeq).reverse_complement():
            if base == "A":
                tale2.rvd += "NI "
            elif base == "T":
                tale2.rvd += "NG "
            elif base == "C":
                tale2.rvd += "HD "
            elif base == "G":
                tale2.rvd += g_RVD

        self.offTargetPairCount = 0

        # Use bitwise operator to compare flag sum to see whether off-target TALEs are on different strands (bad = good cutting ability)
        # or on the same strand (not so bad = FokI domains probably too far apart to cut)
        indivScore = 0

        for (hit1,hit2) in offTargetPairs:
            # Using boolean, count number of offtarget pairs on different strands
            if hit2.flagSum & hit1.flagSum == 0:
                self.diffStrandOffTarget += 1                

            # Using boolean, count number of offtarget pairs on same strand   
            elif hit2.flagSum & hit1.flagSum == 16:               
                self.sameStrandOffTarget += 1
                
            for opt in [hit1.opts, hit2.opts]:
                if opt == "NM:i:0":
                    indivScore += SCORE['INPAIR_OFFTARGET_0']
                if opt == "NM:i:1":
                    indivScore += SCORE['INPAIR_OFFTARGET_1']
                if opt == "NM:i:2":
                    indivScore += SCORE['INPAIR_OFFTARGET_2']
                if opt == "NM:i:3":
                    indivScore += SCORE['INPAIR_OFFTARGET_3']

        # Compute penalties (scores) for off-target hits. Worst = off-target pair, Not so bad = off-target single tale
        self.score = (self.sameStrandOffTarget * SCORE['OFFTARGET_PAIR_SAME_STRAND']) + (self.diffStrandOffTarget * SCORE['OFFTARGET_PAIR_DIFF_STRAND']) + tale1.score + tale2.score + indivScore
        resSites = findRestrictionSites(self.spacerSeq, enzymeCo, minResSiteLen)  
        self.restrictionSites = ";".join(map(lambda x: "%s:%s" % (str(x), ",".join(map(str, resSites[x]))), resSites))
    

    def __str__(self):
        # This creates a tab delimited list of output, with the final column as a semicolon-separated list of REs that cut in the spacer        
        sequence = str(self.tale1.guideSeq) + "*" + self.spacerSeq + "*" + str(self.tale2.guideSeq)

        return "%s\t%s:%s\t%s\t%s\t%s\t%s\t%s\t%s/%s\t%s/%s\t%s/%s\t%s/%s\t%s" % (sequence, self.chrom, self.start, self.tale1.exonNum, self.tale1.rvd, self.tale2.rvd, self.cluster, len(self.offTargetPairs), self.tale1.offTargetsMM[0], self.tale2.offTargetsMM[0], self.tale1.offTargetsMM[1], self.tale2.offTargetsMM[1], self.tale1.offTargetsMM[2], self.tale2.offTargetsMM[2], self.tale1.offTargetsMM[3], self.tale2.offTargetsMM[3], self.restrictionSites)


    def asOffTargetString(self, label, maxOffTargets):
        pairs = []

        # Add any off-target pairs
        if self.offTargetPairs:
            for offTargetPair in self.offTargetPairs:
                pairs.append("%s,%s" % (offTargetPair[0].asOffTargetString(label, maxOffTargets), offTargetPair[1].asOffTargetString(label, maxOffTargets)))
        else:
            pairs.append("")

        pairs = ";".join(pairs)

        return "\n".join([pairs, self.tale1.asOffTargetString("TALE 1", maxOffTargets), self.tale2.asOffTargetString("TALE 2", maxOffTargets)])


class Nickase:
    """ Pair class for 2 Cas9 that are the correct distance apart """
    def __init__(self, tale1, tale2, spacerSeq, spacerSize, offTargetPairs, enzymeCo, maxOffTargets, minResSiteLen):
        self.tale1 = tale1
        self.tale2 = tale2
        self.chrom = tale1.chrom 
        self.strand = tale1.strand
        self.ID = ""
        self.restrictionSites = ""

        # Start of region covered by tale pair
        self.start = tale1.start

        # End of region covered by tale pair
        self.end = tale2.end
        self.spacerSeq = spacerSeq
        self.targetSize = spacerSize
        self.spacerSize = spacerSize
        self.offTargetPairs = offTargetPairs
        self.diffStrandOffTarget = 0
        self.sameStrandOffTarget = 0        


        # Start cluster as -1, but will increment from 1
        self.cluster = -1
        self.spacerStart = tale1.start + tale1.guideSize
        self.spacerEnd = tale2.start - 1

        self.enzymeCo = enzymeCo
        self.strandedGuideSeq = str(self.tale1.guideSeq) + "\n" + self.spacerSeq + "\n" + str(self.tale2.guideSeq)
        self.offTargetPairCount = 0

        # Use bitwise operator to compare flag sum to see whether off-target TALEs are on different strands (bad = good cutting ability)
        # or on the same strand (not so bad = FokI domains probably too far apart to cut)
        indivScore = 0

        for (hit1,hit2) in offTargetPairs:
            # Using boolean, count number of offtarget pairs on different strands
            if hit2.flagSum & hit1.flagSum == 0:
                self.diffStrandOffTarget += 1                
                
            for opt in [hit1.opts, hit2.opts]:
                if opt == "NM:i:0":
                    indivScore += SINGLE_OFFTARGET_SCORE[0]
                if opt == "NM:i:1":
                    indivScore += SINGLE_OFFTARGET_SCORE[1]
                if opt == "NM:i:2":
                    indivScore += SINGLE_OFFTARGET_SCORE[2]
                if opt == "NM:i:3":
                    indivScore += SINGLE_OFFTARGET_SCORE[3]

        # Compute penalties (scores) for off-target hits. Worst = off-target pair, Not so bad = off-target single tale
        self.score = (self.diffStrandOffTarget * SCORE['OFFTARGET_PAIR_DIFF_STRAND']) + tale1.score + tale2.score - indivScore
        resSites = findRestrictionSites(self.spacerSeq, enzymeCo, minResSiteLen)  
        self.restrictionSites = ";".join(map(lambda x: "%s:%s" % (str(x), ",".join(map(str, resSites[x]))), resSites))
    

    def __str__(self):
        # This creates a tab delimited list of output, with the final column as a semicolon-separated list of REs that cut in the spacer        
        sequence = str(self.tale1.guideSeq) + "*" + self.spacerSeq + "*" + str(self.tale2.guideSeq)

        return "%s\t%s:%s\t%s\t%s\t%s\t%s/%s\t%s/%s\t%s/%s\t%s/%s\t%s" % (sequence, self.chrom, self.start, self.tale1.exonNum, self.cluster, len(self.offTargetPairs), self.tale1.offTargetsMM[0], self.tale2.offTargetsMM[0], self.tale1.offTargetsMM[1], self.tale2.offTargetsMM[1], self.tale1.offTargetsMM[2], self.tale2.offTargetsMM[2], self.tale1.offTargetsMM[3], self.tale2.offTargetsMM[3], self.restrictionSites)


    def asOffTargetString(self, label, maxOffTargets):
        pairs = []

        # Add any off-target pairs
        if self.offTargetPairs:
            for offTargetPair in self.offTargetPairs:
                pairs.append("%s,%s" % (offTargetPair[0].asOffTargetString(label, maxOffTargets), offTargetPair[1].asOffTargetString(label, maxOffTargets)))
        else:
            pairs.append("")

        pairs = ";".join(pairs)

        return "\n".join([pairs, self.tale1.asOffTargetString("TALE 1", maxOffTargets), self.tale2.asOffTargetString("TALE 2", maxOffTargets)])

#####################
##
## Functions
##

def scoregRNA(seq, PAM, tail, lookup):
    """ Calculate score from model coefficients. score is 0-1, higher is better """
    score = 0
    if lookup.has_key("Intercept"):
        score = lookup["Intercept"]

    seq = seq[::-1] #we calculate from PAM in a way: 321PAM123

    if lookup.has_key("gc_low"):
        gc = seq[:20].count('G') + seq[:20].count('C')
        if gc < 10:
            score = score + (abs(gc-10) * lookup["gc_low"])
        elif gc > 10:
            score = score + ((gc-10) * lookup["gc_high"])
    
    for i in range(len(seq)):
        key = seq[i] + str(i+1)
        if lookup.has_key(key):
            score += lookup[key]

        if i+1 < len(seq):
            double_key = seq[i] + seq[i+1] + str(i+1)
            if lookup.has_key(double_key):
                score += lookup[double_key]

        if i == 0:
            double_key = PAM[0] + seq[0] + str(0)
            if lookup.has_key(double_key):
                score += lookup[double_key]
    
    for i in range(len(PAM)):
        key = 'PAM' + PAM[i] + str(i+1)
        if lookup.has_key(key):
            score += lookup[key]

        if i+1 < len(PAM):
            double_key = 'PAM' + PAM[i] + PAM[i+1] + str(i+1)
            if lookup.has_key(double_key):
                score += lookup[double_key]
            
    for i in range(len(tail)):
        key = str(i+1) + tail[i]
        if lookup.has_key(key):
            score += lookup[key]

        if i+1 < len(tail):
            double_key = str(i+1) + tail[i] + tail[i+1]
            if lookup.has_key(double_key):
                score += lookup[double_key]
    
    score = 1/(1 + math.e** -score)
    return score


def scoreChari_2015(svmInputFile, svmOutputFile, PAM, genome):
    """ Calculate score from SVM model as in Chari 2015 20-NGG or 20-NNAGAAW, only for hg19 and mm10"""

    model = './models/293T_HiSeq_SP_Nuclease_100_SVM_Model.txt'
    dist = './models/Hg19_RefFlat_Genes_75bp_NoUTRs_SPSites_SVMOutput.txt'
      
    if PAM == 'NGG' and genome == 'mm10':
        model = './models/293T_HiSeq_SP_Nuclease_100_SVM_Model.txt'
        dist = './models/Mm10_RefFlat_Genes_75bp_NoUTRs_SPSites_SVMOutput.txt'
    elif PAM == 'NNAGAAW' and genome == 'hg19':
        model = './models/293T_HiSeq_ST1_Nuclease_100_V2_SVM_Model.txt'
        dist = './models/Hg19_RefFlat_Genes_75bp_NoUTRs_ST1Sites_SVMOutput.txt'
    elif PAM == 'NNAGAAW' and genome == 'mm10':
        model = './models/293T_HiSeq_ST1_Nuclease_100_V2_SVM_Model.txt'
        dist = './models/Mm10_RefFlat_Genes_75bp_NoUTRs_ST1Sites_SVMOutput.txt'

    prog = Popen("./svm_light/svm_classify -v 0 %s %s %s" % (svmInputFile, model, svmOutputFile), shell=True)    
    prog.communicate()

    svmAll = open(dist,'r')
    svmThis = open(svmOutputFile, 'r')

    # first through go all scores and get the max and min
    allData = []
    for line in svmAll:
        line = line.rstrip('\r\n')
        allData.append(float(line))
    svmAll.close()

    scoreArray = []
    for line in svmThis:
        line = line.rstrip('\r\n')
        scoreArray.append(float(line))

    return [ss.percentileofscore(allData, i) for i in scoreArray]


def concatenate_feature_sets(feature_sets):
    '''
    Given a dictionary of sets of features, each in a Pandas.DataFrame,
    concatenate them together to form one big np.array, and get the dimension
    of each set
    Returns: inputs, dim
    Source: Doench 2016
    '''
    assert feature_sets != {}, "no feature sets present"
    F = feature_sets[feature_sets.keys()[0]].shape[0]
    for fset in feature_sets.keys():
        F2 = feature_sets[fset].shape[0]
        assert F == F2, "not same # individuals for features %s and %s" % (feature_sets.keys()[0], fset)

    N = feature_sets[feature_sets.keys()[0]].shape[0]
    inputs = numpy.zeros((N, 0))
    feature_names = []
    dim = {}
    dimsum = 0
    for fset in feature_sets.keys():
        inputs_set = feature_sets[fset].values
        dim[fset] = inputs_set.shape[1]
        dimsum += dim[fset]
        inputs = numpy.hstack((inputs, inputs_set))
        feature_names.extend(feature_sets[fset].columns.tolist())

    return inputs, dim, dimsum, feature_names

def scoreDoench_2016(seq, PAM, tail):
    """ Calculate score from 27-NGG-3 model as in Doench 2016 """
    if len(seq) < 24:
        return 0

    if len(PAM) < 3:
        return 0

    try:
        with open('./models/Doench_2016_model_nopos.pickle', 'rb') as f:
                model= pickle.load(f)

        model, learn_options = model
        learn_options["V"] = 2
        Xdf = pandas.DataFrame(columns=[u'30mer', u'Strand'], data=[[''.join((seq[-24:], PAM, tail[:3])), 'NA']])
        gene_position = pandas.DataFrame(columns=[u'Percent Peptide', u'Amino Acid Cut position'], data=[[-1, -1]])
        feature_sets = feat.featurize_data(Xdf, learn_options, pandas.DataFrame(), gene_position)
        inputs = concatenate_feature_sets(feature_sets)[0]
        return model.predict(inputs)[0]
    except:
        return 0


def gccontent(seq):
        gc = 0

        for i in seq:
                if i =='G' or i == 'g' or i == 'C' or i == 'c':
                        gc = gc + 1
        return( float(gc) / float(len(seq)) )


def get_mismatch_pos(mismatchString):
    mismatches = []

    if mismatchString.isdigit():
        return []

    current = 0
    for c in range(0, len(mismatchString)-1):    

        # If the character is a digit, check if the next character is a digit (>9) and add number to current
        if mismatchString[c].isdigit():
            if mismatchString[c+1].isdigit():
                current += (int(mismatchString[c])*10)
            else:
                current += int(mismatchString[c])
                    
            # if character is a letter, it's a mismatch => add to results
        else:
            mismatches.append(current)
            current += 1

    # Include letters at the end
    if (mismatchString[-1].isalpha()):
        mismatches.append(current)

    return mismatches


def truncateToUTR5(cdsStart, exons, indices):
    """ Truncates the gene to only target 5' UTR """

    endExon = 0

    for exon in range(len(exons)):
        if (cdsStart > exons[exon][1]) and (cdsStart < exons[exon][2]):
            exons[exon][2] = cdsStart
            endExon = exon
            break
        
    return (exons[:(endExon+1)], indices[:(endExon+1)])


def truncateToPROMOTER(strand, exons, indices, bp):
    """ Truncates the gene to only target promoter +-bp TSS """
    
    if (strand == "+"):
        firstExon = exons[0]
        firstExon[2] = firstExon[1] + bp
        firstExon[1] = firstExon[1] - bp
        return ([firstExon], [indices[0]])
    else:
        firstExon = exons[-1]
        firstExon[1] = firstExon[2] - bp
        firstExon[2] = firstExon[2] + bp
        return ([firstExon], [indices[-1]])
    
    return (exons, indices)
    

def truncateToUTR3(cdsEnd, exons, indices):
    """ Truncates the gene to only target 3' UTR """

    startExon = 0

    for exon in range(len(exons)):
        if (cdsEnd > exons[exon][1]) and (cdsEnd < exons[exon][2]):
            exons[exon][1] = cdsEnd
            startExon = exon
        
    return (exons[startExon:], indices[startExon:])


def truncateToSplice(exons, indices):
    """ Truncates the gene to only target splice sites """

    spliceSites = []
    for ind in range(0, len(exons)):
        spliceSites.append([exons[ind][0], exons[ind][1]-1, exons[ind][1]+1])
        spliceSites.append([exons[ind][0], exons[ind][2]-1, exons[ind][2]+1])
        indices.extend([indices[ind],indices[ind]])

    # Remove first and last (i.e. transcription start and termination site)
    return (spliceSites[1:len(spliceSites)-1], indices[1:len(spliceSites)-1])


def truncateToCoding(cdsStart, cdsEnd, exons, indices):
    """ Truncates the gene to only consider the coding region """

    startExon, endExon = 0, len(exons)-1

    # Shortens the coding region to the exons and coordinates between the cds start and cds end
    for exon in range(len(exons)):
        if (cdsStart >= exons[exon][1]) and (cdsStart <= exons[exon][2]):
            exons[exon][1] = cdsStart
            startExon = exon
        
        if (cdsEnd >= exons[exon][1]) and (cdsEnd <= exons[exon][2]):
            # replace the end with the cds end
            exons[exon][2] = cdsEnd
            endExon = exon
    
    if startExon > endExon:
        startExon, endExon = endExon, startExon

    # Shorten list to include exons from cds start to end
    return (exons[startExon:(endExon+1)], indices[startExon:(endExon+1)])

    
def _geneLineToCoord(line, guideSize):
    startSplit = line[1].split(",")
    endSplit = line[2].split(",")
    
    # Remove last (empty) item in list (due to extra comma)
    del startSplit[-1]
    del endSplit[-1]            
    
    startBase = map(int, startSplit)
    endBase = map(int, endSplit)            

    intronSize = [int(startBase[x+1]) - int(endBase[x]) for x in range(len(startBase)-1)]
    intronSize.append(0)
    
    return map(lambda x : [line[0], x[0], x[1], x[2]], zip(startBase, endBase ,intronSize))


def geneToCoord_db(gene, organism, db, guideSize, index):
    """ Gets genomic coordinates for a gene from a database """

    # Try refseq first
    lines = db.execute("SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, refGene r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (organism, gene, gene))

    # Then Ensembl
    if lines == 0: 
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, ensGene r LEFT OUTER JOIN ensemblToGeneName g ON r.name=g.name WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND  (r.name='%s' OR r.name2='%s' OR g.value='%s')" % (organism, gene, gene, gene))

    # Then the general genePred table
    if lines == 0: 
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, r.name, cdsStart, cdsEnd, strand, txStart, txEnd FROM organism o, gpGene r WHERE o.assembly='%s' AND o.organism_id=r.organism_id AND (r.name='%s' OR r.name2='%s')" % (organism, gene, gene))

    # Then wormbase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "ce6" and lines == 0: 
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM sangerGene WHERE (name='%s' OR proteinID='%s')" % (gene, gene))

    # Then flybase. FIX: NO HARDCODED ASSEMBLY!!!
    if organism == "dm3" and lines == 0: 
        lines = db.execute("SELECT chrom, exonStarts, exonEnds, name, cdsStart, cdsEnd, strand, txStart, txEnd FROM flyBaseGene WHERE name='%s'" % (gene))

    if lines == 0:
        sys.stderr.write("The gene name %s was not found in the gene sets for assembly %s. Consider trying an alternative ID (see the instruction page for supported gene identifiers) or using genomic coordinates. If you believe this type of ID should be supported for your organism contact us and we will do our best to support it. \n" % (gene, organism))
        sys.exit(EXIT['GENE_ERROR'])
    
    if lines > 1:
        if index == -1:
            sys.stderr.write("There are multiple transcripts for the gene %s. \n" % gene)
            sys.stderr.write("Please pick a more specific target:\n")
            for line in db:
                sys.stderr.write("%s,%s:%s-%s,%s\n" % (line[3], line[0], line[7], line[8], line[6]))
            sys.exit(EXIT['ISOFORM_ERROR'])
    else:
        index = 1

    for i in range(index):
        line = db.fetchone()

    return (_geneLineToCoord(line, guideSize), line[4], line[5], line[6])


def geneToCoord_file(geneIN, tableFile, guideSize):    
    """ Extracts coordinates of genomic regions to parse for suitable guide binding sites """

    tableR = open(tableFile, 'rb')
    tablereader = csv.DictReader(tableR, delimiter='\t', quoting=csv.QUOTE_NONE)

    # Look in genome table for gene of question
    for row in tablereader:
        if (row['name'] == geneIN or row['name2'] == geneIN):
            line = [row['chrom'], row['exonStarts'], row['exonEnds'], row['cdsStart'], row['cdsEnd'], row['strand']]
            break
    
    tableR.close()
    return (_geneLineToCoord(line[:3], guideSize), line[3], line[4], line[5])

    # Error if gene doesn't exist.
    sys.stderr.write("The gene name %s does not exist in file %s. Please try again.\n" % (geneIN, tableFile))
    sys.exit(EXIT['GENE_ERROR'])

def geneIsoforms(isoform, tableFile):
    
    tableR = open(tableFile, 'rb')
    tablereader = csv.DictReader(tableR, delimiter='\t', quoting=csv.QUOTE_NONE)
    gene_isoforms = set()
    gene = None
    
    for row in tablereader:
        if row['name'] == isoform:
            gene = row['name2']
            break
    
    tableR.close()
    tableR = open(tableFile, 'rb')
    tablereader = csv.DictReader(tableR, delimiter='\t', quoting=csv.QUOTE_NONE)
    for row in tablereader:
        if row['name2'] == gene:
            gene_isoforms.add(row['name'])
    tableR.close()
    
    return gene, isoform, gene_isoforms
    
    
def coordToFasta(regions, fastaFile, outputDir, targetSize, evalAndPrintFunc, indexDir, genome):
    """ Extracts the sequence corresponding to genomic coordinates from a FASTA file """

    sequences = {}
    fastaFile = open(fastaFile, 'w')
    fastaSeq = ""

    for region in regions:
        # Extracts chromosome number and region start and end
        chrom = region[0:region.rfind(':')]
        start = int(region[region.rfind(':')+1:region.rfind('-')])
        finish = int(region[region.rfind('-')+1:])      
        start = max(start, 0)

        # Run twoBitToFa program to get actual dna sequence corresponding to input genomic coordinates
        # Popen runs twoBitToFa program. PIPE pipes stdout.
        prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2> %s/twoBitToFa.err" % (TWOBITTOFA, chrom, start, finish, indexDir, genome, outputDir), stdout = PIPE, shell=True)    

        # Communicate converts stdout to a string
        output = prog.communicate()  
        if (prog.returncode != 0):
            sys.stderr.write("Running twoBitToFa failed\n");
            sys.exit(EXIT['TWOBITTOFA_ERROR']);

        # Ignore other outputs you get from running Popen
        output = output[0]

        # Split string by newlines into a list 
        exons = output.split("\n")        

        # Join together the list without the first line to give just continuous dna sequence of each exon
        dna = ''.join(exons[1:]).upper()

        # Write exon sequences to text file user can open in ApE. exon-intron junctions in lowercase.
        fastaSeq += dna[0].lower()+dna[2:-2]+dna[-1].lower()

        # Add 1 due to BED 0-indexing
        name = "C:%s:%d-%d" % (chrom, start, finish)

        # Loop over exon sequence, write every g-mer into file in which g-mer ends in PAM in fasta format 
        for num in range(0,len(dna)-(targetSize-1)):
            
            if (num - DOWNSTREAM_NUC) > 0:
                start5prim = num - DOWNSTREAM_NUC
            else:
                start5prim = 0
                
            if (num + targetSize + DOWNSTREAM_NUC) > len(dna):
                end3prim = len(dna)
            else:
                end3prim = num + targetSize + DOWNSTREAM_NUC

            downstream_5prim = dna[(start5prim):num] 
            downstream_3prim = dna[(num + targetSize):end3prim]
            if evalAndPrintFunc(name, targetSize, dna[num:(num + targetSize)], num, fastaFile, downstream_5prim, downstream_3prim):
                sequences[name] = dna

    fastaFile.close()
    
    return sequences, fastaSeq


def runBowtie(PAMlength, uniqueMethod_Cong, fastaFile, outputDir, maxOffTargets, indexDir, genome, maxMismatches):
    """ Runs bowtie """

    bowtieResultsFile = '%s/output.sam' % outputDir

    if uniqueMethod_Cong and not ISOFORMS: #When ISOFORMS dna string is not reverse complemented and Cong can't be used
        # the -l alignment mode specifies a seed region to search for the number of mismatches specified with the -n option. Outside of that seed, up to 2 mismatches are searched. 
        # E.g. -l 15 -n 0 will search the first 15 bases with no mismatches, and the rest with up to 3 mismatches
        command = "%s -l %d -n %d -m %d --sam-nohead -k %d %s/%s -f %s -S %s " % (BOWTIE, (PAMlength + 11), maxMismatches, maxOffTargets, maxOffTargets, indexDir, genome, fastaFile, bowtieResultsFile)
    else:
        command = "%s -v %d -m %d --sam-nohead -k %d %s/%s -f %s -S %s " % (BOWTIE, maxMismatches, maxOffTargets, maxOffTargets, indexDir, genome, fastaFile, bowtieResultsFile)
    
    if ISOFORMS: # When ISFORMS we don't check reverse complement 
        command = command + "--norc "
    
    command = command + "2> %s/bowtie.err" % (outputDir)
    
    prog = Popen(command, shell=True)
    prog.wait()

    if (prog.returncode != 0):
        sys.stderr.write("Running bowtie failed\n");
        sys.exit(EXIT['BOWTIE_ERROR']);

    return bowtieResultsFile


def parseBowtie(guideClass, bowtieResultsFile, displayIndices, targets, scoreGC, scoreSelfComp, 
                backbone, replace5prime, maxOffTargets, allowMM, countMM, PAM, scoringMethod=None, 
                genome=None, gene=None, isoform=None, gene_isoforms=None):
    """ Parses bowtie hits and build list of guides"""
    
    currGuide = None
    guideList = []

    with open(bowtieResultsFile, 'rb') as bowtieFile:  
        reader = csv.reader(bowtieFile, delimiter='\t', quoting=csv.QUOTE_NONE)

        for line in reader:
            #  Encountered a new guide RNA (not a new hit for the same guide)
            elements = line[0].split(":") #removes from name 5' and 3' tails
            name = ":".join(elements[0:3])
            if currGuide == None or name != currGuide.name:
                currGuide = guideClass(line[0], line[1], len(line[9]), line[9], scoreGC, scoreSelfComp, 
                                       backbone, PAM, replace5prime, scoringMethod, genome, gene, isoform, gene_isoforms)
                guideList.append(currGuide)

            # Adds hit to off-target list of current guide.
            currGuide.addOffTarget(Hit(line), maxOffTargets, countMM)
        
   
    for guideNum in range(0, len(guideList)):
        for target in range(0, len(targets)):        
            # coordinate for beginning of exon
            targetExon = int(targets[target][ targets[target].rfind(":") + 1:targets[target].rfind("-") ])
            
            # find guides that live in particular exon
            if guideList[guideNum].exonStart == targetExon:                
                # assign them the exon number
                guideList[guideNum].exonNum = displayIndices[target]
        
    return guideList


def parse_primer3_output(target, region, primer3output, primerFastaFile):
    posPattern = re.compile('PRIMER_(\w+)_(\d+)')
    attPattern = re.compile('PRIMER_(\w+)_(\d+)_(\w+)')
    primers = {}
    primerPos = {}
    position, length = 0,0
    
    for line in primer3output.split("\n"):
        if line[0] == "=":
            break 

        label, value = line.split("=")        
        m = attPattern.match(label)
        if m:
            primers[(m.group(2), m.group(1), m.group(3))] = value
        else:
            m = posPattern.match(label)

            if m:
                position, length = value.split(",")            
                
                s, e = int(position), int(position)+int(length)
                if m.group(1) == "RIGHT":
                    s, e = int(position)-int(length)+1, int(position)+1
                primerPos[label] = [s,e]

                primerFastaFile.write(">%s_%s_%s:%s_%s-%s\n%s\n" % (target.ID, m.group(2), m.group(1), region, s, e, primers[(m.group(2), m.group(1), "SEQUENCE")]))

    return primers, primerPos
    

def get_primer_options(options):
    # Parse primer3 options. Update config if known option, otherwise append to primer3 input file
    primerOpt = ""

    if options:
        for opt in options.split(","):
            key, value = opt.split("=")
            if PRIMER3_CONFIG.has_key(key):
                PRIMER3_CONFIG[key] = value
            else:
                primerOpt += opt + "\n"

    return primerOpt

def get_primer_query_sequence_fasta(target, outputDir, flank, fastaSequence):
    s = target.start-flank
    e = target.end+flank
    seqLenBeforeTarget = flank

    if s < 0:        
        seqLenBeforeTarget -= abs(s)
        s = 0

    if e > len(fastaSequence):
        e = len(fastaSequence)

    return fastaSequence[s:e], seqLenBeforeTarget


def get_primer_query_sequence_2bit(target, outputDir, flank, genome, twoBitToFaIndexDir, strand):
    s = target.start-flank
    seqLenBeforeTarget = flank

    if s < 0:        
        seqLenBeforeTarget -= abs(s)
        s = 0

    prog = Popen("%s -seq=%s -start=%d -end=%d %s/%s.2bit stdout 2>> %s/twoBitToFa.err" % (TWOBITTOFA, target.chrom, s, target.end+flank, twoBitToFaIndexDir, genome, outputDir), stdout = PIPE, shell=True)
    output = prog.communicate()  
    
    if (prog.returncode != 0):
        sys.stderr.write("Running twoBitToFa failed\n");
        sys.exit(EXIT['TWOBITTOFA_ERROR']);

    output = output[0].split("\n")
    del(output[0])

    seq = "".join(output)

    return seq, seqLenBeforeTarget


def runBowtiePrimers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets):
    command = "%s -v 0 --best --sam-nohead -k 10 %s/%s -f %s -S %s/primer_results.sam 2> %s/bowtie_primers.err" % (BOWTIE, bowtieIndexDir, genome, primerFastaFileName, outputDir, outputDir)
    prog = Popen(command, shell=True)
    prog.wait()

    if (prog.returncode != 0):
        sys.stderr.write("Running bowtie on primers failed\n");
        sys.exit(EXIT['BOWTIE_PRIMER_ERROR']);

    return parseBowtie(Guide, "%s/primer_results.sam" % outputDir, False, [], [], False, False, None, None, maxOffTargets, None, None, None, None, None)


def make_primers_fasta(targets, outputDir, flanks, genome, limitPrintResults, bowtieIndexDir, fastaSequence, primer3options, guidePadding, enzymeCo, minResSiteLen, geneID, maxOffTargets):
    primers = {}
    primerOpt = get_primer_options(primer3options)

    primerFastaFileName = '%s/primers.fa' % outputDir
    primerFastaFile = open(primerFastaFileName, 'w')
    for i in range(min(limitPrintResults-1, len(targets))):
        target = targets[i]
        seq, seqLenBeforeTarget = get_primer_query_sequence_fasta(target, outputDir, flanks, fastaSequence)
        primer3_output = make_primer_for_target(target, outputDir, seq, seqLenBeforeTarget, primerOpt, guidePadding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start-flanks), min(len(fastaSequence), target.end+flanks))
        target_primers, primerPos = parse_primer3_output(target, region, primer3_output, primerFastaFile)
        primers[target.ID] = target_primers

        # Restriction sites
        restSites = dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen)
        # Sequence for visualization of locus
        dump_locus_sequence(target, outputDir, seq, seqLenBeforeTarget, "+")
        # Genbank file for download
        dump_genbank_file(seq, target, restSites, primerPos, outputDir, geneID, target.start-seqLenBeforeTarget, "+")

    primerFastaFile.close()
       
    primerResults = runBowtiePrimers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets)
    pairPrimers(primers, primerResults, outputDir)


def make_primers_genome(targets, outputDir, flanks, genome, limitPrintResults, bowtieIndexDir, twoBitToFaIndexDir, primer3options, guidePadding, enzymeCo, minResSiteLen, strand, geneID, maxOffTargets):
    primers = {}
    
    primerOpt = get_primer_options(primer3options)

    # RUN PRIMER3 ON TARGET SITES AND CREATE FASTA FILE OF PRIMERS FOR BOWTIE
    primerFastaFileName = '%s/primers.fa' % outputDir
    primerFastaFile = open(primerFastaFileName, 'w')
    for i in range(min(limitPrintResults-1, len(targets))):
        target = targets[i]
        seq, seqLenBeforeTarget = get_primer_query_sequence_2bit(target, outputDir, flanks, genome, twoBitToFaIndexDir, strand)
        primer3_output = make_primer_for_target(target, outputDir, seq, seqLenBeforeTarget, primerOpt, guidePadding)
        region = "%s:%s-%s" % (target.chrom, max(0, target.start-flanks), target.end+flanks)
        target_primers, primerPos = parse_primer3_output(target, region, primer3_output, primerFastaFile)
        primers[target.ID] = target_primers

        # Restriction sites
        restSites = dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen)
        # Sequence for visualization of locus
        dump_locus_sequence(target, outputDir, seq, seqLenBeforeTarget, strand)
        # Genbank file for download
        dump_genbank_file(seq, target, restSites, primerPos, outputDir, geneID, target.start-seqLenBeforeTarget, strand)

    primerFastaFile.close()

    primerResults = runBowtiePrimers(primerFastaFileName, outputDir, genome, bowtieIndexDir, maxOffTargets)
    pairPrimers(primers, primerResults, outputDir)


def dump_restriction_sites(target, seq, flanks, enzymeCo, outputDir, minResSiteLen):
    sites = findRestrictionSites(seq, enzymeCo, minResSiteLen)
    out = [map(lambda x : [str(enzyme), x + target.start-flanks, enzyme.size], sites[enzyme]) for enzyme in sites]
    out = [item for sublist in out for item in sublist]
    out = sorted(out, key=itemgetter(1))

    # Assign tier to avoid overlaps
    siteCount = {}
    tiers = [0] * 23
    for site in out:
        tier = 0

        # count number of sites for each enzyme
        if not siteCount.has_key(site[0]):
            siteCount[site[0]] = 0
        siteCount[site[0]] += 1

        for j in range(len(tiers)):
            if site[1] > tiers[j]:
                tier = j
                tiers[j] = site[1]+site[2]
                break
        site.append(tier)

    # Assign colors depending on uniqueness
    for site in out:
        if siteCount[site[0]] == 1:
            site.append("green")
        else:
            site.append("red")

    outputFile = open("%s/restriction_%s.json" % (outputDir, target.ID), 'w')
    json.dump(out, outputFile)
    outputFile.close()

    return sites


def dump_locus_sequence(target, outputDir, seq, seqLenBeforeTarget, strand):
    if strand == "-":
        seq = str(Seq(seq).complement())
    out = [[target.start-seqLenBeforeTarget, target.end, seq]]
    outputFile = open("%s/locusSeq_%s.json" % (outputDir, target.ID), 'w')
    json.dump(out, outputFile)
    outputFile.close()


def dump_genbank_file(seq, target, restSites, primers, outputDir, geneID, lociStart, strand):
    name= "%s, locus %s" % (geneID, target.ID)
    desc = "CHOPCHOP prediction for gene %s, target %s" % (geneID, target.ID)
    annotation = {"organism" : "Danio rerio", "Target location" : "chrT:1-20"}

    # Genbank file
    genbankFile = open('%s/%s_%s.gb' % (outputDir, geneID, target.ID), 'w')    
    record = SeqRecord(Seq(seq, IUPACAmbiguousDNA()), description=desc, name="CHOPCHOP", id=name)
    record.annotation = annotation

    if target.strand == "+":
        ts = 1
    else:
        ts = -1
        
    record.features.append(SeqFeature(FeatureLocation(target.start-lociStart-1, target.end-lociStart-1, strand=ts), type = "Target"))

    for primer in primers:
        record.features.append(SeqFeature(FeatureLocation(primers[primer][0], primers[primer][1]), type = primer))

    if strand == "-":
        record = record.reverse_complement()

    SeqIO.write(record, genbankFile, "genbank")
    genbankFile.close()

    pass


def pairPrimers(primerAttributes, primerList, outputDir):
    primers = {}

    for primer in primerList:
        guide, primerPairID, side = primer.ID.split("_")

        s = 0
        if side == "RIGHT": s = 1
        if not primers.has_key(guide): primers[guide] = {}
        if not primers[guide].has_key(primerPairID): primers[guide][primerPairID] = [None, None]
        primers[guide][primerPairID][s] = primer

    for guideID in primers:
        guide = primers[guideID]

        att = primerAttributes[int(guideID)]

        outputFile = open("%s/primer_%s.json" % (outputDir, guideID), 'w')
        output = []
        i = 0

        for pairID in guide:
            pair = guide[pairID]

            size = att[(pairID, "PAIR", "PRODUCT_SIZE")]
            ltm = "%.1f" % float(att[(pairID, "LEFT", "TM")])
            rtm = "%.1f" % float(att[(pairID, "RIGHT", "TM")])

            lsq = Seq(att[(pairID, "LEFT", "SEQUENCE")])
            rsq = Seq(att[(pairID, "RIGHT", "SEQUENCE")])

            offTargetPairs = has_Off_targets(pair[0], pair[1], PRIMER_OFF_TARGET_MIN, PRIMER_OFF_TARGET_MIN)
            output.append([ pair[0].chrom, pair[0].start, pair[0].end, pair[1].start, pair[1].end, i, pair[0].strand, "%s" % lsq, "%s" % rsq, len(pair[0].offTargets), len(pair[1].offTargets), len(offTargetPairs), ltm, rtm, size ])

            i += 1

        json.dump(output, outputFile)
        outputFile.close()

 
def make_primer_for_target(guide, outputDir, sequence, seqLenBeforeTarget, primer3options, padding):
    
    template = """PRIMER_SEQUENCE_ID={PRIMER_SEQUENCE_ID:s}
SEQUENCE_TEMPLATE={SEQUENCE_TEMPLATE:s}
SEQUENCE_TARGET={SEQUENCE_TARGET_START:s},{SEQUENCE_TARGET_LEN:s}
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE={PRIMER_OPT_SIZE:s}
PRIMER_MIN_SIZE={PRIMER_MIN_SIZE:s}
PRIMER_MAX_SIZE={PRIMER_MAX_SIZE:s}
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_PRODUCT_SIZE_RANGE={PRODUCT_SIZE_MIN:s}-{PRODUCT_SIZE_MAX:s}
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
"""

    primConfig = PRIMER3_CONFIG.copy()
    primConfig['PRIMER_SEQUENCE_ID'] = str(guide.ID)
    primConfig['SEQUENCE_TEMPLATE'] = sequence
    primConfig['SEQUENCE_TARGET_START'] = str(seqLenBeforeTarget-padding)
    primConfig['SEQUENCE_TARGET_LEN'] = str(guide.targetSize+(2*padding))
    

    primer3InputFile = '%s/%s.primer3Input' % (outputDir, guide.ID)
    f = open(primer3InputFile, 'w')
    f.write(template.format(**primConfig))
    f.write(primer3options)
    f.write("=\n")
    f.close()

    command = "%s < %s 2>> %s/primer3.error" % (PRIMER3, primer3InputFile, outputDir)
    # sys.stderr.write("%s\n" % command)
    prog = Popen(command, stdout = PIPE, shell=True)
    output = prog.communicate()

    if (prog.returncode != 0):
        sys.stderr.write("Running Primer3 failed\n");
        sys.exit(EXIT['PRIMER3_ERROR']);    

    return output[0]


def writeIndividualResults (outputDir, maxOffTargets, sortedOutput, guideSize, mode, totalClusters, limitPrintResults):
    """ Writes each guide and its offtargets into a file """
    
    # Initiate list of lists for each cluster
    clusters = [[] for i in range(totalClusters)]
    
    fileHandler = dict()

    # Limit the number of open files (and results) 
    sortedOutput = sortedOutput[0:min(len(sortedOutput),limitPrintResults-1)] 

    for i in range(len(sortedOutput)):
        current = sortedOutput[i]
        current.ID = i+1

        # Create new file if not already opened
        if not current.ID in fileHandler:
            resultsFile = '%s/%s.offtargets' % (outputDir, current.ID)            
            fileHandler[current.ID] = open(resultsFile, 'w')
        f = fileHandler[current.ID]

        # Add the current TALE pair to the appropriate list in the list of lists, depending on its cluster number            
        if mode == TALENS or mode == NICKASE:
            clusterID = current.cluster
            clusters[clusterID-1].append(current)             
            
        offTargets = current.asOffTargetString("", maxOffTargets)
        if not offTargets:
            offTargets = "There are no predicted off-targets."

        f.write(str(current.strandedGuideSeq)+"\n"+offTargets+"\n")
                
    for clust in clusters:
        if len(clust) == 0:
            continue
        bestInCluster = clust[0]

        for member in clust[1:]:
            # Write the other cluster members to file
            fileHandler[bestInCluster.ID].write("%s*%s*%s,%s:%s,%s,%s,%s/%s,%s/%s,%s/%s,%s/%s;" % (member.tale1.guideSeq, member.spacerSeq, member.tale2.guideSeq, member.chrom, member.start, member.tale1.exonNum, len(member.offTargetPairs), member.tale1.offTargetsMM[0], member.tale2.offTargetsMM[0], member.tale1.offTargetsMM[1], member.tale2.offTargetsMM[1], member.tale1.offTargetsMM[2], member.tale2.offTargetsMM[2], member.tale1.offTargetsMM[3], member.tale2.offTargetsMM[3]))

        fileHandler[bestInCluster.ID].write("\n"+current.restrictionSites+"\n")

    for fh in fileHandler.values():
        fh.close()
        
    return clusters            


def findRestrictionSites(sequence, enzymeCompany, minSize=1):
    # Take spacerSeq as DNA input for restriction site search
    mySeq = Seq(sequence, IUPACAmbiguousDNA())

    # Restricts enzyme possibilities to NEB enzymes. Can ultimately change to any supplier.
    rb = RestrictionBatch(first=[], suppliers=[enzymeCompany])        

    # Filter binding sites shorter than given length
    rb = filter(lambda x: len(x) > minSize, rb)

    # Determine which restriction enzymes cut in the sequence provided
    analyze = Analysis(rb, mySeq)
    return analyze.with_sites()


def comaprePAM(basePAM, baseDNA):
    if basePAM == "N":
        return True

    if basePAM == baseDNA:
        return True

    if basePAM == "W" and (baseDNA == "A" or baseDNA == "T"):
        return True

    if basePAM == "S" and (baseDNA == "C" or baseDNA == "G"):
        return True

    if basePAM == "M" and (baseDNA == "A" or baseDNA == "C"):
        return True

    if basePAM == "K" and (baseDNA == "G" or baseDNA == "T"):
        return True

    if basePAM == "R" and (baseDNA == "A" or baseDNA == "G"):
        return True

    if basePAM == "Y" and (baseDNA == "C" or baseDNA == "T"):
        return True

    if basePAM == "B" and baseDNA != "A":
        return True

    if basePAM == "D" and baseDNA != "C":
        return True

    if basePAM == "H" and baseDNA != "G":
        return True

    if basePAM == "V" and baseDNA != "T":
        return True

    return False

#####################
##
## CPF1 SPECIFIC FUNCTIONS
##


def eval_CPF1_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM):
    """ Evaluates an k-mer as a potential Cpf1 target site """

    gLen = guideSize-len(PAM)
    revCompPAM = str(Seq(PAM).reverse_complement())
    dna = Seq(dna)
    
    add = True
    for pos in range(len(PAM)):
        if comaprePAM(PAM[pos], dna[pos]):
            continue 
        else:
            add = False
            break
                
    if add:
        dna = dna.reverse_complement() if not ISOFORMS else dna
        fastaFile.write('>%s_%d-%d:%s:%s:+\n%s\n' % (name, num, num+guideSize, downstream5prim, downstream3prim, dna))
        return True
    
    add = True
    
    for pos in range(len(PAM)):
        if comaprePAM(revCompPAM[pos], dna[gLen + pos]):
            continue 
        else:
            add = False
            break
    
    if add and not ISOFORMS:
        #on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
        fastaFile.write('>%s_%d-%d:%s:%s:-\n%s\n' % (name, num, num+guideSize, Seq(downstream3prim).reverse_complement(), Seq(downstream5prim).reverse_complement(), dna))
        return True

    return False


#####################
##
## CRISPR SPECIFIC FUNCTIONS
##


def eval_CRISPR_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed, PAM):
    """ Evaluates an k-mer as a potential CRISPR target site """

    gLen = guideSize-len(PAM)
    revCompPAM = str(Seq(PAM).reverse_complement())
    dna = Seq(dna)

    if str(dna[0:2]) in allowed:
        add = True
        for pos in range(len(PAM)):
            if comaprePAM(PAM[pos], dna[gLen + pos]):
                continue 
            else:
                add = False
                break
            
        # in order to control the number of mismatches to search in the last 8 or 3 bps, 
        # need to reverse complement so the seed region can be at the start
        # rather than end of the sequence
        # not in isoforms case as we don't search reverse complement
        if add:
            dna = dna.reverse_complement() if not ISOFORMS else dna
            fastaFile.write('>%s_%d-%d:%s:%s:+\n%s\n' % (name, num, num+guideSize, downstream5prim, downstream3prim, dna))
            return True

    if str(dna[-2:].reverse_complement()) in allowed and not ISOFORMS:
        add = True

        for pos in range(len(PAM)):
            if comaprePAM(revCompPAM[pos], dna[pos]):
                continue 
            else:
                add = False
                break

        if add:
            #on the reverse strand seq of 5' downstream becomes 3' downstream and vice versa
            fastaFile.write('>%s_%d-%d:%s:%s:-\n%s\n' % (name, num, num+guideSize, Seq(downstream3prim).reverse_complement(), Seq(downstream5prim).reverse_complement(), dna))
            return True

    return False


def sort_CRISPR_guides(guides):
    """ Sort pairs according to score  """
    return sorted(guides, key=attrgetter('score'))


#####################
##
## TALEN SPECIFIC FUNCTIONS
##


def pairTalens(taleList, fastaSeq, guideSize, taleMinDistance, taleMaxDistance, enzymeCo, maxOffTargets, g_RVD, minResSiteLen):
    pairs = []

    for i in range(len(taleList)-1):
        tale1 = taleList[i]

        # FIX: Only start looking for pair when > 30 - 36 spacer+length of i-TALE (modified for 17-mers and 18-mers)
        for j in range(i+1, len(taleList)-1):
            tale2 = taleList[j]

            # This will finish the search for more pairs if we are out of range
            if tale1.start + taleMaxDistance < tale2.start:
                break 

            # If pairs are within the correct range of one another and one tale begins with T, while the other ends with A...
            # elif tale1.start + taleMinDistance < tale2.start and tale1.guideSeq[0] == "T" and tale2.guideSeq[guideSize-1] == "A" and tale1.strand == "+" and tale2.strand == "-":
            elif tale1.start + taleMinDistance < tale2.start and tale1.guideSeq[0] == "T" and tale2.guideSeq[guideSize-1] == "A":

                # EDV: Are all these find calls faster than a regular expression?
                pos = tale1.name.find('_')
                exon1 = tale1.name[:pos]
                exonSeq = fastaSeq[exon1]

                # Make sure the two TALENs are on the same "slice", only a problem for overlapping padding regions
                pos2 = tale2.name.find('_')
                exon2 = tale2.name[:pos2]
                if exon1 != exon2:
                    continue

                # The coordinates of the tale within the exon e.g. 128-143
                tale1coords = tale1.name[pos+1:]

                # Just the second coordinate, corresponding to the end of the first tale e.g. 143
                tale1End = int(tale1coords[tale1coords.find('-')+1:])

                # The coordinates of the tale within the exon e.g. 160-175
                tale2coords = tale2.name[tale2.name.find('_')+1:]                            

                # Just the first coordinate, corresponding to the beginning of the second tale e.g. 160
                tale2Start = int(tale2coords[:tale2coords.find('-')])

                # sequence of spacer between end of tale1 and beginning of tale2
                spacerSeq = exonSeq[tale1End:tale2Start]

                spacerSize = len(spacerSeq)                         

                # if spacerSize < 3:
                #      sys.stderr.write("(%s)  (%s)\n" % (tale1.name, tale2.name))
                #      sys.stderr.write("(%s)  (%s)\n" % (e1, e2))
                #      sys.stderr.write("%s-%s\n" % (tale1End, tale2Start))
                #      sys.stderr.write("%s\t%s\t%s\n" % (tale1.guideSeq, spacerSeq, tale2.guideSeq))
                #      sys.exit()

                # Calculates off-target pairs for tale1 and tale2 (see below)
                offTargetPairs = has_Off_targets(tale1, tale2, TALEN_OFF_TARGET_MIN, TALEN_OFF_TARGET_MAX)

                # Makes tale1 and tale2 into a Pair object, and adds to list of Pair objects
                pairs.append(Pair(tale1, tale2, spacerSeq, spacerSize, offTargetPairs, enzymeCo, maxOffTargets, g_RVD, minResSiteLen))

    return pairs


def pairCas9(taleList, fastaSeq, guideSize, taleMinDistance, taleMaxDistance, enzymeCo, maxOffTargets, minResSiteLen, offtargetMaxDist):
    pairs = []

    for i in range(len(taleList)-1):
        tale1 = taleList[i]

        # FIX: Only start looking for pair when > 30 - 36 spacer+length of i-TALE (modified for 17-mers and 18-mers)
        for j in range(i+1, len(taleList)-1):
            tale2 = taleList[j]

            if tale1.start + taleMaxDistance < tale2.start:
                continue 

            elif tale1.start + taleMinDistance < tale2.start and tale1.strand != tale2.strand:

                # EDV: Are all these find calls faster than a regular expression?
                pos = tale1.name.find('_')
                exon1 = tale1.name[:pos]
                exonSeq = fastaSeq[exon1]

                # Make sure the two TALENs are on the same "slice", only a problem for overlapping padding regions
                pos2 = tale2.name.find('_')
                exon2 = tale2.name[:pos2]
                if exon1 != exon2:
                    continue

                # The coordinates of the tale within the exon e.g. 128-143
                tale1coords = tale1.name[pos+1:]

                # Just the second coordinate, corresponding to the end of the first tale e.g. 143
                tale1End = int(tale1coords[tale1coords.find('-')+1:])

                # The coordinates of the tale within the exon e.g. 160-175
                tale2coords = tale2.name[tale2.name.find('_')+1:]                            

                # Just the first coordinate, corresponding to the beginning of the second tale e.g. 160
                tale2Start = int(tale2coords[:tale2coords.find('-')])

                # sequence of spacer between end of tale1 and beginning of tale2
                spacerSeq = exonSeq[tale1End:tale2Start]

                spacerSize = len(spacerSeq)                         

                # Calculates off-target pairs for tale1 and tale2 (see below)
                offTargetPairs = has_Off_targets(tale1, tale2, taleMinDistance-guideSize, offtargetMaxDist)

                # Makes tale1 and tale2 into a Pair object, and adds to list of Pair objects
                pairs.append(Nickase(tale1, tale2, spacerSeq, spacerSize, offTargetPairs, enzymeCo, maxOffTargets, minResSiteLen))

    return pairs


def has_Off_targets(tale1, tale2, offTargetMin, offTargetMax):
    """ Returns the number of off-targets for a pair of TALENs (10-24bp apart) """

    offTargetPairs = []  

    # Calls sort function to sort off-targets by chromosome and chromosome position. Bowtie ranks them according to quality of hit
    tale1.sort_offTargets()
    tale2.sort_offTargets()
  
    ### FIX: Eivind to write this code properly. Include a way to step backwards, so as not to miss any hits. Need to make a queue..?
    for i in range(len(tale1.offTargets)):
        hit1 = tale1.offTargets[i]

        for j in range(len(tale2.offTargets)):
            hit2 = tale2.offTargets[j]          

            # Determines whether 2 tales are on the same chromosome and 10-24 bp apart.
            if hit2.chrom == hit1.chrom and offTargetMin <= abs(hit2.start-hit1.start) <= offTargetMax:
                offTargetPairs.append([hit1, hit2])

    return offTargetPairs


def clusterPairs(pairs):
    """ Clusters paired sequences according to overlap, so user knows which TALE pairs are redundant """

    # Sets the starting pair of TALEs to be compared to
    first = pairs[0]
    cluster = 1
    first.cluster = cluster
    inCluster = 0

    # Compares each TALE pair to previous pair in list to see whether redundant. Assigns cluster number accordingly
    for i in range(1,len(pairs)):
        cur = pairs[i]      
        prev = pairs[i-1]
        
        # Specifically, compares location of spacer (by comparing location of tales) to see whether there is overlap, and therefore TALE pairs are redundant
        if ((cur.spacerStart <= prev.spacerEnd) and (cur.spacerEnd >= prev.spacerStart) and inCluster < PRIMER_OFF_TARGET_MIN):

            # Checks whether spacer overlap is >50%, in which case tale pairs are somewhat redundant, and should be part of the same cluster
            # if (cur.spacerEnd - prev.spacerStart+1) > (0.9*(cur.spacerSize+prev.spacerSize)):
                # Adds current cluster to pairs list
            cur.cluster = cluster
            inCluster += 1
        else:
            # If not redundant, increase cluster number
            cluster += 1
            cur.cluster = cluster
            inCluster = 0
    
    return (cluster, pairs)


def eval_TALENS_sequence(name, targetSize, dna, num, fastaFile, downstream5prim, downstream3prim):
    """ Evaluates an N-mer as a potential TALENs target site """
    del downstream5prim, downstream3prim
    found = False

    # One TALE must start with a T and the other must end in an A, such that each on its respective strand begins with a T
    if dna[0] == "T":
        # dna = Seq(dna).reverse_complement()    
        fastaFile.write('>%s_%d-%d\n%s\n' % (name, num, num+targetSize, dna))
        found = True
    elif dna[-1] == "A":
        fastaFile.write('>%s_%d-%d\n%s\n' % (name, num, num+targetSize, dna))
        found = True

    return found


def sort_TALEN_pairs(pairs):
    """ Sort pairs according to score and cluster """

    return sorted(pairs, key=attrgetter('score', 'cluster'))


#####################
##
## JASON visualization
##


def complement(sequence):
    return sequence.translate(string.maketrans("ACGT", "TGCA"))


def FastaToViscoords(sequences, strand):
    """ Makes the exons in 'sequences' array generated in coordToFasta json readable for visualization"""
    exonstart = []
    exonend = []
    exonsequence = []

    for exon in sequences:
        # sys.stderr.write("%s\n" % exon)
        exonlist = exon.split(':')
        exoncoord = exonlist[2].split('-')
        exonstart.append(exoncoord[0])
        exonend.append(exoncoord[1])
        seq = sequences[exon]
        if strand == "-":
            seq = complement(seq)

        exonsequence.append(seq)

    return zip(exonstart, exonend, exonsequence)

#####################
##
## MAIN
##

def getAllowedFivePrime(allowed):

    newAllowed = []
    for el in allowed.split(","):
        if el[0] == 'N' and el[1] == 'N':
            return ("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
        elif el[0] == 'N':
            newAllowed.extend(["A"+el[1], "C"+el[1], "G"+el[1], "T"+el[1]])
        elif el[1] == 'N':
            newAllowed.extend([el[0]+"A", el[0]+"C", el[0]+"G", el[0]+"T"])
        else:
            newAllowed.append(el)

    return dict(zip(newAllowed, [True] * len(newAllowed)))


def coordsToJson(coords, cdsStart, cdsEnd, strand):
    newCoords = []

    for coord in coords:
        if coord[2] <= cdsStart:
            coord.append(VIS_WIDTH['UTR'])
            coord.append(strand)
            newCoords.append(coord)
        elif coord[1] >= cdsEnd: 
            coord.append(VIS_WIDTH['UTR'])
            coord.append(strand)
            newCoords.append(coord)
        else:
            if cdsStart > coord[1] and cdsStart < coord[2]:
                newCoords.append([coord[0], coord[1], cdsStart, 0, VIS_WIDTH['UTR'], strand])
#                newCoords.append([coord[0], cdsStart, coord[2], coord[3], VIS_WIDTH['CDS']])
                coord[1] = cdsStart


            if cdsEnd > coord[1] and cdsEnd < coord[2]:
#                newCoords.append([coord[0], coord[1], cdsEnd, 0, VIS_WIDTH['CDS']])
                newCoords.append([coord[0], cdsEnd, coord[2], coord[3], VIS_WIDTH['UTR'], strand])
                coord[2] = cdsEnd
                coord[3] = 0

            newCoords.append([coord[0], coord[1], coord[2], coord[3], VIS_WIDTH['CDS'], strand])

#           if coord[1] >= cdsStart and coord[2] <= cdsEnd:
#                newCoords.append([coord[0], coord[1], coord[2], coord[3], VIS_WIDTH['CDS']])


#    sys.stderr.write("NEWCOORDS: %s\n" % newCoords)
    return newCoords


def parseTargets(targetString, genome, use_db, data, padSize, targetRegion, exonSubset, promoter_bp):
    targets = []
    displayIndices = []
    visCoords = []
    targetSize = 0
    target_strand = None

    pattern = re.compile("(([\.\w]+):)?([\.\,\d]+)\-([\.\,\d]+)")
    isCoordinate = pattern.match(targetString)

    if isCoordinate: 
        
        if ISOFORMS:
                sys.stderr.write("--isoforms is not yet ready with coordinate search.\n")
                sys.exit(EXIT['ISOFORMS_ERROR'])
                
        # Target specified as coordinate
        if target_strand == None:
            target_strand = "+"
        elif target_strand == "-":
            sys.stderr.write("All targets must be on the same strand.\n")
            sys.exit(EXIT['GENE_ERROR'])


        ## CHROMOSOME LOCATION
        chrom = isCoordinate.group(2)
        i = 1

        for target in targetString.split(";"):
            m = pattern.match(target)

            if m:
                if m.group(2) != None and chrom != m.group(2):
                    sys.stderr.write("Can't target regions on separate chromosomes (%s != %s).\n" % (chrom, m.group(1)))
                    sys.exit(EXIT['GENE_ERROR'])

                startPos = m.group(3)
                endPos = m.group(4)

                startPos = startPos.replace(",","").replace(".","")
                endPos = endPos.replace(",","").replace(".","")

                startPos = int(startPos)
                endPos = int(endPos)
                targetSize += endPos - startPos + 1

                if (startPos >= endPos):
                    sys.stderr.write("Start position (%s) must be smaller than end position (%s)\n" % (startPos, endPos))
                    sys.exit(EXIT['GENE_ERROR'])
                    
                targets.append("%s:%s-%s" % (chrom, max(0, startPos-padSize), endPos+padSize))
                visCoords.append([chrom, startPos, endPos, 0, 10, "+"])
                displayIndices.append(i)
                i += 1
            else:
                sys.stderr.write("Unknown format: %s\n" % (target))
                sys.exit(EXIT['GENE_ERROR'])
                
    else:
        # Check whether index is provided. For cases where you have multiple entries with same name.
        index = -1
        pattern = re.compile("(.+):::(\d+)")
        m = pattern.match(targetString)
        if m:
            targetString = m.group(1)
            index = int(m.group(2))

        ## GENE / TRANSCRIPT
        if use_db:
            (visCoords, cdsStart, cdsEnd, strand) = geneToCoord_db(targetString, genome, data, padSize, index)
            if ISOFORMS:
                sys.stderr.write("--isoforms is not yet ready with database search.\n")
                sys.exit(EXIT['ISOFORMS_ERROR'])
        else:
            (visCoords, cdsStart, cdsEnd, strand) = geneToCoord_file(targetString, data, padSize)

        if target_strand == None:
            target_strand = strand
        elif target_strand != strand:
            sys.stderr.write("All targets must be on the same strand.\n")
            sys.exit(EXIT['GENE_ERROR'])

        # Subset on exons
        coords = copy.deepcopy(visCoords)
        visCoords = coordsToJson(visCoords, cdsStart, cdsEnd, strand)

        if strand == "-":
            coords.reverse()

        coords, displayIndices = subsetExons(exonSubset, coords)

        if strand == "-":
            coords.reverse()
            displayIndices.reverse()

        # Truncate to region
        if targetRegion == "CODING":
            coords, displayIndices = truncateToCoding(cdsStart, cdsEnd, coords, displayIndices)
        elif targetRegion == "UTR5":
            if strand == "+":
                coords, displayIndices = truncateToUTR5(cdsStart, coords, displayIndices)
            else:
                coords, displayIndices = truncateToUTR3(cdsEnd, coords, displayIndices)
        elif targetRegion == "PROMOTER":
                coords, displayIndices = truncateToPROMOTER(strand, coords, displayIndices, promoter_bp)
        elif targetRegion == "UTR3":
            if strand == "+":
                coords, displayIndices = truncateToUTR3(cdsEnd, coords, displayIndices)
            else:
                coords, displayIndices = truncateToUTR5(cdsStart, coords, displayIndices)
        elif targetRegion == "SPLICE":
            coords, displayIndices = truncateToSplice(coords, displayIndices)
        elif targetRegion != "WHOLE":
            sys.stderr.write("Unknown region: %s\n" % targetRegion);
            sys.exit(EXIT['PYTHON_ERROR']);

        for exon in coords:
            targetSize += exon[2] - exon[1] + 1

        # Pad since can bind outside exons
        targets.extend(map(lambda x : "%s:%s-%s" % (x[0], x[1]-padSize, x[2]+padSize), coords))

    if targetSize > TARGET_MAX:
        sys.stderr.write("Search region is too large (%s nt). Maximum search region is %s nt.\n" % (targetSize, TARGET_MAX))
        sys.exit(EXIT['GENE_ERROR'])

    return (targets, displayIndices, visCoords, target_strand)


def parseFastaTarget(fastaFile, candidateFastaFile, targetSize, evalAndPrintFunc):
    """ Parse a FASTA file as input for targeting """

    seq_name = "sequence"
    sequence = ""
    fastaFile = open(fastaFile, 'r')
    for line in fastaFile:
        if line[0] == ">":
            seq_name = '_'.join(line[1:].strip().split())
        else:
            sequence += line.rstrip()        
    fastaFile.close()

    
    name = "%s:0-%s" % (seq_name, len(sequence))
    idName = "C:" + name 
    sequence = sequence.upper()
    sequence = "".join(sequence.split())

    dna_pattern = re.compile(r'([^ACGTNacgtn])')    
    if dna_pattern.search(sequence):
        sys.stderr.write("Input sequence contains illegal characters.\n")
        sys.exit(EXIT['GENE_ERROR']);

        
    sequences = {}
    candidateFastaFile = open(candidateFastaFile, 'w')

    # Loop over sequence, write every k-mer into file in which k-mer ends in as PAM in fasta format 
    for num in range(0,len(sequence)-(targetSize-1)):
        
        if (num - DOWNSTREAM_NUC) > 0:
                start5prim = num - DOWNSTREAM_NUC
        else:
                start5prim = 0
                
        if (num + targetSize + DOWNSTREAM_NUC) > len(sequence):
                end3prim = len(sequence)
        else:
                end3prim = num + targetSize + DOWNSTREAM_NUC

        downstream_5prim = sequence[(start5prim):num] 
        downstream_3prim = sequence[(num + targetSize):end3prim]
        
        if evalAndPrintFunc(idName, targetSize, sequence[num:(num + targetSize)], num, candidateFastaFile, downstream_5prim, downstream_3prim):
            sequences[idName] = sequence

    return (sequences, [name], [1], [[seq_name, 1, len(sequence), 0, 20, "+"]], sequence, "+")


def hyphen_range(s):
    """ Takes a range in form of "a-b" and generate a list of numbers between a and b inclusive.
    Also accepts comma separated ranges like "a-b,c-d,f" will build a list which will include
    Numbers from a to b, a to d and f"""

    s = "".join(s.split()) #removes white space
    r = set()

    for x in s.split(','):
        t = x.split('-')
        if len(t) not in [1,2]: raise SyntaxError("Range is not properly formatted: "+ s)
        if len(t)==1: 
            r.add(int(t[0])) 
        else: 
            r.update(set(range(int(t[0]),int(t[1])+1)))

    l = list(r)
    l.sort()

    return l


def subsetExons(exons, targets):
    indices = range(len(targets)+1)[1:]

    if exons:
        indices = hyphen_range(exons)

        for index in indices:
            if int(index) > len(targets):
                sys.stderr.write("That exon does not exist\n")
                sys.exit(EXIT['PYTHON_ERROR'])
                
        targets = [targets[int(i)-1] for i in indices]      # indices is a list of exon numbers -1 e.g. exon 2 is [1]

    return (targets, indices)


def connect_db(database_string):
    import MySQLdb
    
    m = re.compile("(\w+):(\w+)@(\w+)/(\w+)").match(database_string)

    if not m:
        sys.stderr.write("Wrong syntax for connection string: %s\n" % database_string)
        sys.exit(EXIT["DB_ERROR"])

    try:
        db = MySQLdb.connect(user=m.group(1), passwd=m.group(2), host=m.group(3), db=m.group(4))
    except:
        sys.stderr.write("Could not connect to database\n")
        sys.exit(EXIT['DB_ERROR'])

    return db


def getMismatchVectors(pam, gLength, cong):
    
    allowed = [True] * (gLength -len(pam))
    count = [True] * (gLength -len(pam))

    if cong:
        allowed = [True] * 9 + [False] * (gLength -len(pam) -9)

    for char in pam:
        count.append(False)
        if char == "N":
            allowed.append(True)
        else:
            allowed.append(False)

    return allowed, count


def getCpf1MismatchVectors(pam, gLength):
    
    allowed = [True] * (gLength -len(pam))
    count = [True] * (gLength -len(pam))

    for char in pam[::-1]:
        count.insert(0, False)
        if char == "N":
            allowed.insert(0,True)
        else:
            allowed.insert(0,False)

    return allowed, count


def mode_select(var, index, MODE):
    """ Selects a default depending on mode for options that have not been set """

    if var is not None:
        return var
        
    if MODE == CRISPR:
        return CRISPR_DEFAULT[index]

    elif MODE == TALENS:
        return TALEN_DEFAULT[index]
    
    elif MODE == CPF1:
        return CPF1_DEFAULT[index]

    elif MODE == NICKASE:
        return NICKASE_DEFAULT[index]

    sys.stderr.write("Unknown model %s\n" % MODE)
    sys.exit(EXIT['PYTHON_ERROR'])


def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", help="Target genes or regions", metavar="TARGET_REGIONS")
    parser.add_argument("-r", "--gRVD", default="NH ", dest="g_RVD", action="store_const", const="NN ",  help="Use RVD 'NN' instead of 'NH' for guanine nucleotides. 'NH' appears to be more specific than 'NN' but the choice depends on assembly kit.")
    parser.add_argument("-D", "--database", help="Connect to a chopchop database to retrieve gene: user_name:passwd@host/database", metavar="DATABASE", dest="database")   
    parser.add_argument("-e", "--exon", help="Comma separated list of exon indices. Only find sites in this subset. ", metavar="EXON_NUMBER", dest="exons")   
    parser.add_argument("-TP", "--targetPromoter", default=200, type=int, help="how many bp to target upstream and downstream of TSS")
    parser.add_argument("-G", "--genome", default="danRer7", metavar="GENOME", help="The genome to search.") 
    parser.add_argument("-g", "--guideSize", default=None, type=int, metavar="GUIDE_SIZE", help="The size of the guide RNA.")
    parser.add_argument("-c", "--scoreGC", default=None, action="store_false", help="Score GC content. True for CRISPR, False for TALENs.")
    parser.add_argument("-SC", "--noScoreSelfComp", default=None, action="store_false", help="Do not penalize self-complementarity of CRISPR.")
    parser.add_argument("-BB", "--backbone", default=None, metavar="BACKBONE", help="Penalize self-complementarity versus backbone regions (comma-separated list, same strand as guide). Requires -C.")
    parser.add_argument("-R5", "--replace5P", default=None, metavar="REPLACE_5P", help="Replace bases from 5' end (with e.g. 'GG') ")  ## FIX: AT THE MOMENT THIS IS ONLY APPLIES TO FOLDING/SELF-COMPL
    parser.add_argument("-t", "--target", default="CODING", dest="targetRegion", help="Target the whole gene CODING/WHOLE/UTR5/UTR3/SPLICE. Default is CODING.")
    parser.add_argument("-T", "--MODE", default=1, type=int, choices=[1, 2, 3, 4], help="Set mode (int): default is Cas9 = 1, Talen = 2, Cpf1 = 3, Nickase = 4")
    parser.add_argument("--taleMin", default=14, type=int, help="Minimum distance between TALENs. Default is 14.")  # 14 + 18(length of TALE) = 32
    parser.add_argument("--taleMax", default=20, type=int, help="Maximum distance between TALENs. Default is 20.")  # 20 + 18(length of TALE) = 38
    parser.add_argument("--nickaseMin", default=10, type=int, help="Minimum distance between TALENs. Default is 10.") 
    parser.add_argument("--nickaseMax", default=31, type=int, help="Maximum distance between TALENs. Default is 31.")
    parser.add_argument("--offtargetMaxDist", default=100, type=int, help="Maximum distance between offtargets for Nickase. Default is 100.")
    parser.add_argument("-f", "--fivePrimeEnd", default="NN", metavar="FIVE_PRIME_END", help="Specifies the requirement of the two nucleotides 5' end of the CRISPR guide: A/C/G/T/N. Default: NN.")
    parser.add_argument("-n", "--enzymeCo", default="N", metavar="ENZYME_CO", help="The restriction enzyme company for TALEN spacer.")
    parser.add_argument("-R", "--minResSiteLen", type=int, default=4, help="The minimum length of the restriction enzyme.")
    parser.add_argument("-v", "--maxMismatches", default=3, type=int, choices=[0, 1, 2, 3], metavar="MAX_MISMATCHES", help="The number of mismatches to check across the sequence.")
    parser.add_argument("-m", "--maxOffTargets", metavar="MAX_HITS", help="The maximum number of off targets allowed.")
    parser.add_argument("-M", "--PAM", type=str, help="The PAM motif.")
    parser.add_argument("-o", "--outputDir", default="./", metavar="OUTPUT_DIR", help="The output directory. Default is the current directory.")
    parser.add_argument("-F", "--fasta", default=False, action="store_true", help="Use FASTA file as input rather than gene or genomic region.")
    parser.add_argument("-p", "--padSize", default=-1, type=int, help="Extra bases searched outside the exon. Defaults to the size of the guide RNA for CRISPR and TALEN + maximum spacer for TALENS.")
    parser.add_argument("-P", "--makePrimers", default=False, action="store_true", help="Designes primers using Primer3 to detect mutation.")
    parser.add_argument("-3", "--primer3options", default=None, help="Options for Primer3. E.g. 'KEY1=VALUE1;KEY2=VALUE2'")
    parser.add_argument("-A", "--primerFlanks", default=300, type=int, help="Size of flanking regions to search for primers.")
    parser.add_argument("-a", "--guidePadding", default=20, type=int, help="Minimum distance of primer to target site.")
    parser.add_argument("-O", "--limitPrintResults", default=1000, dest="limitPrintResults", help="The number of results to print extended information for. Default 1000.")
    parser.add_argument("-w", "--uniqueMethod_Cong", default=False, dest="uniqueMethod_Cong", action="store_true", help="A method to determine how unique the site is in the genome: allows 0 mismatches in last 15 bp.")
    parser.add_argument("-J", "--jsonVisualize", default=False, action="store_true", help="Create files for visualization with json.")
    parser.add_argument("--scoringMethod", default="G_20", type = str, choices=["XU_2015", "DOENCH_2014", "DOENCH_2016", "MORENO_MATEOS_2015", "CHARI_2015", "G_20", "ALL"], help="Scoring used for Cas9 and Nickase. Default is G_20")
    parser.add_argument("--rm1perfOff", default = False, action="store_true", help="For fasta input, don't score one off-target without mismatches.")
    parser.add_argument("--isoforms", default = False, action="store_true", help="Search for offtargets on the transcriptome.")
    args = parser.parse_args()
    
    # set isoforms to global as it is influencing many steps
    global ISOFORMS
    ISOFORMS = args.isoforms

    # Add TALEN length
    args.taleMin += 18
    args.taleMax += 18

    # Set mode specific parameters if not set by user
    args.scoreGC = mode_select(args.scoreGC, "SCORE_GC", args.MODE)
    args.scoreSelfComp = mode_select(args.noScoreSelfComp, "SCORE_FOLDING", args.MODE)
    args.PAM = mode_select(args.PAM, "PAM", args.MODE)
    args.guideSize = mode_select(args.guideSize, "GUIDE_SIZE", args.MODE) + len(args.PAM)
    args.maxMismatches = mode_select(args.maxMismatches, "MAX_MISMATCHES", args.MODE)
    args.maxOffTargets = mode_select(args.maxOffTargets, "MAX_OFFTARGETS", args.MODE)

    # Add TALEN length
    args.nickaseMin += args.guideSize
    args.nickaseMax += args.guideSize

    if args.scoreSelfComp:
        if args.backbone:
            tmp = args.backbone.strip().split(",")
            args.backbone = [str(Seq(el).reverse_complement()) for el in tmp]
        else:
            args.backbone = []

    # Pad each exon equal to guidesize unless
    if args.padSize != -1:
        padSize = args.padSize
    else:
        if args.MODE == TALENS:
            padSize = args.taleMax
        elif args.MODE == NICKASE:
            padSize = args.nickaseMax
        elif args.MODE == CRISPR or args.MODE == CPF1:
            padSize = args.guideSize

    # Set default functions for different modes  
    if args.MODE == CRISPR:
        # Set mismatch checking policy
        (allowedMM, countMM) = getMismatchVectors(args.PAM, args.guideSize, args.uniqueMethod_Cong)
        allowed = getAllowedFivePrime(args.fivePrimeEnd)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CRISPR_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed=allowed, PAM=args.PAM)
        guideClass = Cas9 if not ISOFORMS else Guide
        sortOutput = sort_CRISPR_guides
    elif args.MODE == CPF1:
        (allowedMM, countMM) = getCpf1MismatchVectors(args.PAM, args.guideSize)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CPF1_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, PAM=args.PAM)
        guideClass = Guide
        sortOutput = sort_CRISPR_guides
    elif args.MODE == TALENS:
        (allowedMM, countMM) = getMismatchVectors(args.PAM, args.guideSize, None)
        guideClass = Guide
        evalSequence = eval_TALENS_sequence
        sortOutput = sort_TALEN_pairs
    elif args.MODE == NICKASE:
        (allowedMM, countMM) = getMismatchVectors(args.PAM, args.guideSize, args.uniqueMethod_Cong)
        allowed = getAllowedFivePrime(args.fivePrimeEnd)
        evalSequence = lambda name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim: eval_CRISPR_sequence(name, guideSize, dna, num, fastaFile, downstream5prim, downstream3prim, allowed=allowed, PAM=args.PAM)
        guideClass = Cas9
        sortOutput = sort_TALEN_pairs

    # Connect to database if requested
    if args.database:
        cdb = connect_db(args.database)
        db = cdb.cursor()
        use_db = True
    else:
        db = "%s/%s.gene_table" % (GENE_TABLE_INDEX_DIR, args.genome)
        use_db = False

    ## Create output directory if it doesn't exist
    if not os.path.isdir(args.outputDir):
        os.mkdir(args.outputDir)
    
    candidateFastaFile = '%s/sequence.fa' % args.outputDir
    if args.fasta:
        sequences, targets, displayIndices, visCoords, fastaSequence, strand = parseFastaTarget(args.targets, candidateFastaFile, args.guideSize, evalSequence)        
    else:
        targets, displayIndices, visCoords, strand = parseTargets(args.targets, args.genome, use_db, db, padSize, args.targetRegion, args.exons, args.targetPromoter)
        sequences, fastaSequence = coordToFasta(targets, candidateFastaFile, args.outputDir, args.guideSize, evalSequence, TWOBIT_INDEX_DIR, args.genome)
    
    gene, isoform, gene_isoforms = geneIsoforms(args.targets, db) if ISOFORMS else (None, None, set())
    
    ## Converts genomic coordinates to fasta file of all possible k-mers
    if len(sequences) == 0:
        sys.stderr.write("No target sites\n")
        sys.exit()
        
    # Run bowtie and get results
    bowtieResultsFile = runBowtie(len(args.PAM), args.uniqueMethod_Cong, candidateFastaFile, args.outputDir, 
                                  int(args.maxOffTargets), ISOFORMS_INDEX_DIR if ISOFORMS else BOWTIE_INDEX_DIR, 
                                  args.genome, int(args.maxMismatches))
    results = parseBowtie(guideClass, bowtieResultsFile, displayIndices, targets, args.scoreGC, args.scoreSelfComp, 
                          args.backbone, args.replace5P, args.maxOffTargets, allowedMM, countMM, args.PAM, 
                          args.scoringMethod, args.genome, gene, isoform, gene_isoforms)  # TALENS: MAKE_PAIRS + CLUSTER

    if args.rm1perfOff and args.fasta:
        for guide in results:
            if guide.offTargetsMM[0] > 0:
                guide.score = guide.score - SINGLE_OFFTARGET_SCORE[0]
                
    if ISOFORMS:
        for guide in results:
            other_isoforms = guide.gene_isoforms.remove(guide.isoform) if guide.isoform in guide.gene_isoforms else guide.gene_isoforms
            other_off_isoforms = guide.offTargetsIso[0].remove(guide.isoform) if guide.isoform in guide.offTargetsIso[0] else guide.offTargetsIso[0]
            guide.conserved = other_isoforms == other_off_isoforms

    if (args.scoringMethod == "CHARI_2015" or args.scoringMethod == "ALL") and (args.PAM == "NGG" or args.PAM == "NNAGAAW") and (args.genome == "hg19" or args.genome == "mm10") and not ISOFORMS:
        try:
            #make file to score
            svmInputFile = '%s/chari_score.SVMInput.txt' % args.outputDir
            svmOutputFile = '%s/chari_score.SVMOutput.txt' % args.outputDir
            encoding = defaultdict(str)
            encoding['A'] = '0001'
            encoding['C'] = '0010'
            encoding['T'] = '0100'
            encoding['G'] = '1000'

            svmFile = open(svmInputFile, 'w')

            for guide in results:
                seq = guide.downstream5prim + guide.strandedGuideSeq[:-len(guide.PAM)]
                PAM = guide.strandedGuideSeq[-len(guide.PAM):]
                sequence = (seq[-20:] + PAM).upper()
                x = 0
                tw = '-1'
                # end index
                if len(sequence)==27:
                    endIndex = 22
                else:
                    endIndex = 21
                
                while x < endIndex:
                    y = 0
                    while y < 4:
                        tw = tw + ' ' + str(x+1) + str(y+1) + ':' + encoding[sequence[x]][y]
                        y += 1
                    x += 1
                svmFile.write(tw + '\n')
            svmFile.close()
            newScores = scoreChari_2015(svmInputFile, svmOutputFile, args.PAM, args.genome)

            for i, guide in enumerate(results):
                guide.CoefficientsScore["CHARI_2015"] = newScores[i]
                if args.scoringMethod == "CHARI_2015":
                    guide.score = guide.score - (guide.CoefficientsScore["CHARI_2015"]/100) * SCORE['COEFFICIENTS']
        except:
            pass


    if args.MODE == CRISPR or args.MODE == CPF1 or ISOFORMS:
        cluster = 0
    elif args.MODE == TALENS:
        pairs = pairTalens(results, sequences, args.guideSize, int(args.taleMin), int(args.taleMax), args.enzymeCo, args.maxOffTargets, args.g_RVD, args.minResSiteLen)

        if (not len(pairs)):
            sys.stderr.write("No TALEN pairs could be generated for this region.\n")
            sys.exit(EXIT['GENE_ERROR'])

        if args.rm1perfOff and args.fasta:
            for pair in pairs:
                if pair.diffStrandOffTarget > 0:
                    pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]
                if pair.sameStrandOffTarget > 0:
                    pair.score = pair.score - SCORE["OFFTARGET_PAIR_SAME_STRAND"]

        cluster, results = clusterPairs(pairs)
        
    elif args.MODE == NICKASE:
        pairs = pairCas9(results, sequences, args.guideSize, int(args.nickaseMin), int(args.nickaseMax), args.enzymeCo, args.maxOffTargets, args.minResSiteLen, args.offtargetMaxDist)

        if (not len(pairs)):
            sys.stderr.write("No Cas9 nickase pairs could be generated for this region.\n")
            sys.exit(EXIT['GENE_ERROR'])

        if args.rm1perfOff and args.fasta:
            for pair in pairs:
                if pair.diffStrandOffTarget > 0:
                    pair.score = pair.score - SCORE["OFFTARGET_PAIR_DIFF_STRAND"]

        cluster, results = clusterPairs(pairs)

    # Sorts pairs according to score/penalty and cluster 
    if strand == "-" and not ISOFORMS:
        results.reverse()

    sortedOutput = sortOutput(results)

    # Write individual results to file
    listOfClusters = writeIndividualResults(args.outputDir, args.maxOffTargets, sortedOutput, args.guideSize, args.MODE, cluster, args.limitPrintResults)

    if args.makePrimers:
        if args.fasta:
            make_primers_fasta(sortedOutput, args.outputDir, args.primerFlanks, args.genome, args.limitPrintResults, BOWTIE_INDEX_DIR, fastaSequence, args.primer3options, args.guidePadding, args.enzymeCo, args.minResSiteLen, "sequence", args.maxOffTargets)
        else:
            make_primers_genome(sortedOutput, args.outputDir, args.primerFlanks, args.genome, args.limitPrintResults, BOWTIE_INDEX_DIR, TWOBIT_INDEX_DIR, args.primer3options, args.guidePadding, args.enzymeCo, args.minResSiteLen, strand, args.targets, args.maxOffTargets)


    ## Print results 
    resultCoords = []
    
    if ISOFORMS:
        print "Rank\tTarget sequence\tGenomic location\tExon\tGene\tIsoform\tGC content (%)\tSelf-complementarity\tMM0\tMM1\tMM2\tMM3\tConserved\tIsoformsMM0\tIsoformsMM1\tIsoformsMM2\tIsoformsMM3"
        for i in range(len(sortedOutput)):
            print "%s\t%s" % (i+1, sortedOutput[i])            
            resultCoords.append([sortedOutput[i].start, sortedOutput[i].score, sortedOutput[i].guideSize, sortedOutput[i].strand])
    else:                                                                         
        if args.MODE == CRISPR:
            common_header = "Rank\tTarget sequence\tGenomic location\tExon\tStrand\tGC content (%)\tSelf-complementarity\tMM0\tMM1\tMM2\tMM3"
            if args.scoringMethod == "ALL":
                print(common_header + "\tXU_2015\tDOENCH_2014\tDOENCH_2016\tMORENO_MATEOS_2015\tCHARI_2015\tG_20")
            else:
                print(common_header + "\tEfficiency")
            for i in range(len(sortedOutput)):
                print "%s\t%s" % (i+1, sortedOutput[i])            
                resultCoords.append([sortedOutput[i].start, sortedOutput[i].score, sortedOutput[i].guideSize, sortedOutput[i].strand])
                
        elif args.MODE == CPF1:
            print "Rank\tTarget sequence\tGenomic location\tExon\tStrand\tGC content (%)\tSelf-complementarity\tMM0\tMM1\tMM2\tMM3"
            for i in range(len(sortedOutput)):
                print "%s\t%s" % (i+1, sortedOutput[i])            
                resultCoords.append([sortedOutput[i].start, sortedOutput[i].score, sortedOutput[i].guideSize, sortedOutput[i].strand])
    
        elif args.MODE == TALENS or args.MODE == NICKASE:
            
            if args.MODE == TALENS:
                print "Rank\tTarget sequence\tGenomic location\tExon\tTALE 1\tTALE 2\tCluster\tOff-target pairs\tOff-targets MM0\tOff-targets MM1\tOff-targets MM2\tOff-targets MM3\tRestriction sites\tBest ID"
            else:
                print "Rank\tTarget sequence\tGenomic location\tExon\tCluster\tOff-target pairs\tOff-targets MM0\tOff-targets MM1\tOff-targets MM2\tOff-targets MM3\tRestriction sites\tBest ID"
            finalOutput = []
            for cluster in listOfClusters:  ## FIX: WHY ARE THERE EMPTY CLUSTERS???
                if len(cluster) == 0:
                    continue
    
                finalOutput.append(cluster[0])
    
            sortedFinalOutput = sortOutput(finalOutput)
            resultCoords = [[j+1, sortedFinalOutput[j].spacerStart, sortedFinalOutput[j].score, sortedFinalOutput[j].spacerSize, sortedFinalOutput[j].strand, sortedFinalOutput[j].ID, sortedFinalOutput[j].tale1.start, sortedFinalOutput[j].tale2.end] for j in range(len(sortedFinalOutput))]
    
            for i in range(len(sortedFinalOutput)):
                print "%s\t%s\t%s" % (i+1,sortedFinalOutput[i], sortedFinalOutput[i].ID)
    
    # Print gene annotation files

    # FASTA file
    geneFile = open('%s/gene_file.fa' % args.outputDir, 'w')
    geneFile.write(">%s\n" % args.targets)
    geneFile.write(fastaSequence)
    geneFile.close()


    # Visualize with json
    if args.jsonVisualize:
        strand = visCoords[0][5] # set based on first exon

        # Coordinates for gene
        visCoordsFile = open('%s/viscoords.json' % args.outputDir, 'w')
        visCoords = sorted(visCoords,  key=itemgetter(1))
        json.dump(visCoords, visCoordsFile)

        # Coordinates for sequence
        seqvis = FastaToViscoords(sequences, strand)
        seqvisFile = open('%s/seqviscoords.json' % args.outputDir, 'w')
        json.dump(seqvis, seqvisFile)
     
        # Coordinates for cutters
        cutCoord_file = open('%s/cutcoords.json' % args.outputDir, 'w') 

        cutcoords = []
        for i in range(len(resultCoords)):
            el = []

            if args.MODE == CRISPR or args.MODE == CPF1:                
                el.append(i+1)
                el.extend(resultCoords[i])
            elif args.MODE == TALENS or args.MODE == NICKASE:
                el.extend(resultCoords[i])
                
            cutcoords.append(el)


        # Put bars at different heights to avoid overlap
        tiers = [0] * 23
        sortedCoords = sorted(cutcoords, key=itemgetter(1))
        for coord in sortedCoords:

            t = 0
            for j in range(len(tiers)):
                if coord[1] > tiers[j]:
                    t = j
                    tiers[j] = coord[1]+coord[3]
                    break

            coord.append(t)
    
        json.dump(cutcoords, cutCoord_file)

        info = open("%s/run.info" % args.outputDir, 'w')
        info.write("%s\t%s\t%s\t%s\t%s\n" % ("".join(args.targets), args.genome, args.MODE, args.uniqueMethod_Cong, args.guideSize))
        info.close()

if __name__ == '__main__':
    main()

