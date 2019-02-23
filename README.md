# CHOPCHOP script
#### There exists website version of this tool with specifically designed visualization of guides under [chopchop.cbu.uib.no](http://chopchop.cbu.uib.no)

#### This repository is open sourced as specified in the LICENSE file. It is Apache License 2.0.

#### About:
CHOPCHOP is a python script that allows quick and customizable design of guide RNA. We support selecting target sites for CRISPR/Cas9, CRISPR/Cpf1 or TALEN with wide range of customization. We even support Cas13 for isoform targeting.


#### Prerequisites:
Please, create separate environment for CHOPCHOP using `virtualenv` to easly manage dependencies.  

- [Python](https://www.python.org/download/) - We operate on 2.7
- [Biopython module](http://biopython.org/wiki/Download "Biopython module download")
- Python libraries: pandas, numpy, pickle, scipy, argparse, MySQLdb (if you intend to use our SQL database)
- Additionally Python library [scikit-learn==0.18.1](https://pypi.python.org/pypi/scikit-learn/0.18.1#downloads) if you want to make use of ```--scoringMethod DOENCH_2016```, otherwise latest version is ok. This older version is required because models from Doench et al. 2016 have been saved with this particular version.
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
    * Fill name with extension ".gene_table' e.g. danRer10.gene_table
    * Get output
- [Download](http://hgdownload.soe.ucsc.edu/downloads.html) *.2bit compressed genome
    * Select organism in complete annotation sets section
    * Select Full data set
    * download *.2bit file
- Create fasta version of genome by running twoBitToFa on *.2bit file
- [Make bowtie compressed version of genome](http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer) using your new *.fasta file
- Copy config.json and create config_local.json file, replace paths with your own for 
*.2bit genome files, bowtie (*.ewbt) genome files and *.gene_table files
- Make sure all these files and programs have proper access rights
- Have fun using CHOPCHOP as a script  

When using CHOPCHOP for targeting isoforms (e.g. Cas13a, option ```--isoforms```) one needs also to create 
transcript version of .fasta file, where each new  description line is describing name of the isoform and 
following sequence is the sequence of the isoform, reverse complemented if necessary. This can be easily 
achieved with bedtools getfasta e.g.  
```bedtools getfasta -fi danRer10.fa -bed danRer10.bed -name -fo danRer10.transciptome.fa -s -split```  
Bowtie indexes of transcriptome files should also be created. In this situation all possible guides 
for given isoform will be created, and mismatches will be checked against the transcriptome. Additionally 
column "Constitutive" will indicate True when guide is conserved in the whole family of isoforms of the gene, 
and False otherwise. Also, column "IsoformsMM0" will contain names of the isoforms (of target isoform gene 
family) that are also targeted by the guide with 0 mismatches. Column MM0 will contain number of off-targets 
with 0 mismatches, but without counting of-targets on the same isoform family.

It is also possible for ```--isoforms``` mode to use base pairing probability as efficiency, for this you need to 
set up a folder with .mt files. To create .mt files you can use ViennaRNA package and run:
  ```
 RNAplfold < your_tx_file.fa
  ```
Afterwards use their `mountain.pl` script to creat .mt files, for e.g.
```bash
#!/bin/bash  
for filename in $1/*.ps; do  
	./mountain.pl < $filename > "$2/$(basename "$filename" _dp.ps).mt"  
done  
```
And finally set he path to those .mt files in config_local.json. It is important to name folder the same way you named 
the genome files, and transcript files. You can check how its set up on chopchop website or download files from there 
directly.
  
 - All available chopchop genomes are [downloadable](http://chopchop.cbu.uib.no/bin/genomes/)
 - Latest SQL database is also available in the folder above, named e.g. chopchop_dev_20180427.sql, you can use 
 this databse instead of .gene_table files
 - Isoform folder with isoform related indexes and .mt files is [here](http://chopchop.cbu.uib.no/bin/genomes/isoforms/)
  

#### Run example:
List gRNAs using default values for CRIPR/Cas9 for `chr10:1000000-1001000`, with genome named danRer10 and put results in directory temp:
  
  ```
  ./chopchop.py -G danRer10 -o temp -Target chr10:1000000-1001000
  ```

List gRNAs using default values for CRIPR/Cas9 for gene NM_144906, with genome named hg19 and put results in directory temp:
  
  ```
  ./chopchop.py -G hg19 -o temp -Target NM_144906
  ```

To imitate output from the website you can use template below to get similar results, but closely examine options on the
website and adjust settings below. Many options on the website do not have exactly the same name, but read the descriptions
of each parameter. For the website [valenvm.cbu.uib.no](valenvm.cbu.uib.no) it is possible to get almost all parameters of your query 
(without database) by simply visiting "http://valenvm.cbu.uib.no/results/YOUR_RUN_ID/query.json". 

  ```
  chopchop.py -J -P -T 1 -M NGG --maxMismatches 3 -g 20  --scoringMethod DOENCH_2016  -f NN --backbone AGGCTAGTCCGT --replace5P GG  -G hg38  -t CODING -n N  -R 4  -3 'PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60' -A 290 -a 20 --rm1perfOff -o temporary/ your_gene_name > temporary/results.txt 2> temporary/python.err
  ```
  
Also, to make sure that chopchop will prefer to create output in case where one of the machine learning alghoritms fails 
(we are trying to use other authors code without many changes) we are outputing zeros in the Efficiency columns when 
that happens. Therefore if you see Efficiency column being all zeros, it probably means you might not have properly 
installed packages required by that alghoritm e.g. DOENCH_2016 or KIM_2018.  

In cases where you see differences in off-targets between website and your own installation of CHOPCHOP, make sure that 
you are using the same verson of the genome and bowtie indexes as the website. 
Genomes and indexes are available [here](http://chopchop.cbu.uib.no/bin/genomes/).  

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
  
  Design Cas13 guides for selected transcripts of tb gene in Zebrafish:
  ```
  ./chopchop_query.py --gene_names ENSDART00000007204.8,ENSDART00000160271.1,ENSDART00000157768.1 -G danRer10 -g 27 -M H -T 3 --isoforms -o /tb/ENSDARG00000039806
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