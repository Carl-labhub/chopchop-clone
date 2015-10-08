# CHOPCHOP

#### Prerequisites:
- [Python](https://www.python.org/download/ "Python latest version")
- [Biopython module](http://biopython.org/wiki/Download "Biopython module download")
- [Bowtie](http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.1/ "Bowtie download")
- [twoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/ "twoBitToFa download")
- [primer3](http://primer3.sourceforge.net/releases.php "primer3 download") - optional
- Add bowtie and twoBitToFa (primer3 also if you want primers designed) to PATH by editing etc/profile
- CHOPCHOP will need a [table](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) to look up genomic cooridinates
    * Select organism and assembly 
    * Select group: Genes and Gene Predictions
    * Select track: RefSeq Genes 
    * Select table: refFlat 
    * Select region: genome
    * Select output format: all fields from selected table
    * Fill name with extension ".UCSC_table' eg. danRer7.UCSC_table
    * Get output
- [Download](http://hgdownload.soe.ucsc.edu/downloads.html) *.2bit compressed genome
    * Select organism in complete annotation sets section
    * Select Full data set
    * download *.2bit file
- Create fasta version of genome by running twoBitToFa on *.2bit file
- [Make bowtie compressed version of genome](http://bowtie-bio.sourceforge.net/manual.shtml#the-bowtie-build-indexer) using your new *.fasta file
- Create a genomes directory in chopchop directory and put there: *.2bit genome file, bowtie (*.ewbt) genome files and *.UCSC_table file
- Make sure all these files and programs have proper acces rights

#### Run examples:
- Calculate gRNAs using default values for CRIPR/Cas9 for chr10:1000000-1001000, with genome named danRer10 and put results in directory temp:
  
  ./chopchop.py -G danRer10 -o temp chr10:1000000-1001000
