# CHOPCHOP script
#### There exists website version of this tool with specifically designed visualization of guides under [chopchop.cbu.uib.no](http://chopchop.cbu.uib.no)

#### This repository is open sourced as specified in the LICENSE file. It is Apache License 2.0.

#### About:
CHOPCHOP is a python script that allows quick and customizable design of guide RNA. We support selecting target sites for CRISPR/Cas9, CRISPR/Cpf1 or TALEN with wide range of customization. We even support C2c2 for isoform targeting.


#### Prerequisites:
- [Python](https://www.python.org/download/) - We operate on 2.7
- [Biopython module](http://biopython.org/wiki/Download "Biopython module download")
- Python libraries: pandas, numpy, pickle, scipy
- Additionally Python library [skcit-learn==0.18.1](https://pypi.python.org/pypi/scikit-learn/0.18.1#downloads) if you want to make use of ```--scoringMethod DOENCH_2016```, otherwise latest version is ok. This older version is required because models from Doench et al. 2016 have been saved with this particular version.
- Additionally Python library keras with theano backend if you want to use KIM et al 2018 model for Cpf1 efficiency
- [Bowtie](http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.1/ "Bowtie download") - included, but may require compilation for your operating system
- [twoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/ "twoBitToFa download") - included
- [svm_light](http://svmlight.joachims.org/ "svm_light download") - included, but may require compilation for your operating system, necessary only with option ```--scoringMethod CHARI_2015```, only working for mm10 and hg19 genomes with NGG or NNAGAAW PAMs, otherwise returns zeros
- [primer3](http://primer3.sourceforge.net/releases.php "primer3 download") - included
- CHOPCHOP script will need a [table](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) to look up genomic coordinates if you want to supply names of the genes rather than coordinates. To get example genePred table:
    * Select organism and assembly 
    * Select group: Genes and Gene Predictions
    * Select track: RefSeq Genes or Ensemble Genes 
    * Select table: refFlat or ensGene
    * Select region: genome
    * Select output format: all fields from selected table
    * Fill name with extension ".gene_table' eg. danRer10.gene_table
    * Get output
- [Download](http://hgdownload.soe.ucsc.edu/downloads.html) *.2bit compressed genome
    * Select organism in complete annotation sets section
    * Select Full data set
    * download *.2bit file
- Create fasta version of genome by running twoBitToFa on *.2bit file
- [Make bowtie compressed version of genome](http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer) using your new *.fasta file
- In script chopchop.py (lines 43-45) set paths to your *.2bit genome files, bowtie (*.ewbt) genome files and *.gene_table files
- Make sure all these files and programs have proper access rights
- Have fun using CHOPCHOP as a script  

When using CHOPCHOP for targeting isoforms (eg. C2c2, option ```--isoforms```) one needs also to create transcript version of .fasta file, where each new  description line is describing name of the isoform and following sequence is the sequence of the isoform, reverse complemented if necessary. This can be easily achieved with bedtools getfasta eg. ```bedtools getfasta -fi danRer10.fa -bed danRer10.bed -name -fo danRer10.transciptome.fa -s -split```  
Bowtie indexes of transcriptome files should also be created. In this situation all possible guides for given isoform will be created, and mismatches will be checked against the transcriptome. Additionally column "Conserved" will indicate True when guide is conserved in the whole family of isoforms of the gene, and False otherwise. Also, column "IsoformsMM0" will contain names of the isoforms (of target isoform gene family) that are also targeted by the guide with 0 mismatches. Column MM0 will contain number of off-targets with 0 mismatches, but without counting of-targets on the same isoform family.


#### Run example:
List gRNAs using default values for CRIPR/Cas9 for chr10:1000000-1001000, with genome named danRer10 and put results in directory temp:
  
  ```
  ./chopchop.py -G danRer10 -o temp chr10:1000000-1001000
  ```

List gRNAs using default values for CRIPR/Cas9 for gene NM_144906, with genome named hg19 and put results in directory temp:
  
  ```
  ./chopchop.py -G hg19 -o temp NM_144906
  ```

#### Additionally we include:  
```control_guides.py``` - script to find CRSIPR Cas9/Cpf1 guides that do not map to selected genome, follow specific GC content, have no self-complementarity or complementarity to the backbone and are filtered for supplied list of restriction sites  

  ```
  ./control_guides.py /path/to/bowtie/index/of/the/genome --PAM TTN --type Cas9 --how_many 400 --g_len 20 --restrict BbsI,EcoRI
  ```

```chopchop_query.py``` - script to find guides using CHOPCHOP for a list of genes or all genes in the genePred file. You can use all chopchop.py script options normally, they will be passed along in each query. 
  

  For two selected genes:  
  ```
  ./chopchop_query.py --gene_names NM_144906,NM_001265875 -G hg19 -o temp2
  ```

  For all genes in genePred table, with Cpf1 option:  
  ```
  ./chopchop_query.py --genePred /full/path/to/genPred/hg19.gene_table -G hg19 -o temp -T 3
  ```
  
  Design C2c2 guides for selected transcripts of tb gene in Zebrafish:
  ```
  ./chopchop_query.py --gene_names ENSDART00000007204.8,ENSDART00000160271.1,ENSDART00000157768.1 -G danRer10 -g 27 -M C --isoforms -o /tb/ENSDARG00000039806
  ```  
  
#### Explore different options of our CHOPCHOP scripts:
  ```
  ./chopchop.py --help
  ```  

  ```
  ./control_guides.py --help
  ```  

  ```
  ./chopchop_query.py --help
  ```