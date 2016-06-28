## VCF Converter to FASTA format with some statistics
Next Generation Sequencing (NGS) generates a huge amount of data. One of the formats that allows the storage of information about sequences and their annotations is Variant Call Format (VCF). This thesis describes the main features of VCF file format, in particular focusing on the description of information about polymorphic sites in sequences. They have an important role in population studies, among other things. The key element of this dissertation is a converter from VCF to multi-FASTA format. It is a script written in Python, which can be executed from command line. Therefore it is possible to execute the program in different pipelines. The appended manual includes the description of all available options as well as the explanation of the converter’s functionalities.

### General information
The converter is written in Python 2.7.6. <br />
It is possible to analyze VCF files for single individuals, but also for an entire population. <br />
The VCF file has to be in a valid format, e.g. for single individual : <br />
> \#\#fileformat=VCFv4.1 <br />
> \#\#FILTER=<ID=LowQual,Description="Low quality"> <br />
> \#\#FORMAT=<ID=DP,Number=1,Type=Integer> <br />
> \#\#INFO=<ID=AN,Number=1,Type=Integer> <br />
> \#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1549 <br />
> ace 1 . C . 76.23 . AN=2;DP=21 GT:DP 0|0:21 <br />
> ace 2 . T . 76.23 . AN=2;DP=21 GT:DP 0|0:21 <br />
> ace 3 . A . 82.23 . AN=2;DP=23 GT:DP 0|0:23 <br />
> ace 4 . C T 55.77 . AN=2;DP=23 GT:DP 0/1:22 <br />
> ace 5 . C . 82.23 . AN=2;DP=23 GT:DP 0|0:23 <br />
> ace 6 . C . 79.23 . AN=2;DP=23 GT:DP 0|0:23 <br />

For the correct execution of the script, the header line starting with `#CHROM` has to be present. `1549` in the example is the id of my individual (its presence is important too).

By default, if the difference in DP value between two neighbour nucleotides is greater than 80%, the letter `n` is transmitted to the output file. This value can be modified in the script according to your required precision.

### Usage
The only thing you need is to download the vcf_converter.py and parser.py files.  
Then run:  
`python vcf_converter.py [options]`

For more information please run:  
`python vcf_converter.py –h`  
This will output the following:
```
Usage: vcf_converter.py [options]  
This is a converter from VCF format to FASTA.  
Options:  
-h, --help                            show this help message and exit  
-t INPUT_TYPE, --type=INPUT_TYPE      chose 'folder' or 'file' to be executed  
-f FILENAME, --file=FILENAME          read data from FILENAME  
-i INPUT_PATH, --input=INPUT_PATH     read input path  
-o OUTPUT_PATH, --output=OUTPUT_PATH  read output path  
```

####Options
- **-t INPUT_TYPE, --type=INPUT_TYPE**: it is compulsory to choose one of the available options. Type `-t file` if you want to process a single file; type `-t folder` if you want to process an entire folder of files; in this case only files in vcf format are analyzed. If inside the folder there are other kind of files or other folders (even containing other vcf files) they will be ignored.
- **-f FILENAME, --file=FILENAME**: this option is necessary only if you chose the *file* option as INPUT_TYPE; example of valid FILENAME: *abh_1545.vcf*
- **-i INPUT_PATH, --input=INPUT_PATH**: in this option the relative path is taken into consideration (the path corresponding to `pwd` in Linux). It **can't** have the backslash at the beginning. This option can be skipped if you want to process files in the same directory of the main scripts (this path is set by default).
- **-o OUTPUT_PATH, --output=OUTPUT_PATH**: destination path for the output files. The rules are the same as for the INPUT_PATH.

### Output files
- FASTA file, containing the sequences of the two alleles (A and B) for each individual for each gene; if a gene hasn't been already analyzed, a new file called *gene\_name.fasta* will be created. The alleles' sequences of further individuals will be appended to existing files.  
Example of a file for two analyzed individuals (file name: *abh.fasta*):
```
>abh_1622_A
CAnGAAGTAGACTAGCTACAGGCTCAATTTGATAAGAGGAGGTCAACTCGGTCTCCTTGTAGGGCAGTGGAACTCTCCTTTATACCCAGAGAGTAAATGGGCACCAGTATCACATAAGGGAGATGAGCTCTTCTATGGGTACAGTGTTTTCATGGGGTAGCTGTTTGTGGTCTCATCACAGAGCATTGGTCTTTTATGTTGTTTTGAGGATACCCAGGTAGCACAGTTGTTGAGTCTAGTGTCTAAAGGAGCAGAAACGTTATATGGTCTGTGAACCAGTATCATTCAGATGATGATGTCTTCTAGATGTTGAGCTCTGCATGTATGTAAACCCAGCCTTTGAGTAATGCCTCTTCCTTAGGCTAGTTTTGGTTTTACCCTAGAAGTTGGAAGGGTGCGTAGATATTTTTTTAACTGAAATTTGTTCTTAAATACAAAAAGATAACATAGAATTAAGCTGATCTTTAGGACAGATAACATTGATGCAAAACTGTGAAGTA
>abh_1622_B
CAnGAAGTAGACTAGCTACAGGCTCAATTTGATAAGAGGAGGTCAACTCGGTCTCCTTGTAGGGCAGTGGAACTCTCCTTTATACCCAGAGAGTAAATGGGCACCAGTATCACATAAGGGAGATGAGCTCTTCTATGGGTACAGTGTTTTCATGGGGTAGCTGTTTGTGGTCTCATCACAGAGCATTGGTCTTTTATGTTGTTTTGAGGATACCCAGGTAGCACAGTTGTTGAGTCTAGTGTCTAAAGGAGCAGAAACGTTATATGGTCTATGAACCAGTATCATTCAGATGATGATGTCTTCTAGATGTTGAGCTCTGCATGTATGTAAACCCAGCCTTTGAGTAATACCTCTTCCTTAGGCTAGTTTTGGTTCTACCCTAGAAGTTGGAAGGGTGCGTAGATATTTTTTTAACTGAAATTTGTTCTTAAATACAAAAAGATAACATAGAATTAAGCTGATCTTTAGGACAGATAACATTGATGCAAAACTGTGAAGTA
>abh_1627_A
CAnGAAGTAGACTAGCTACAGGCTCAATTTGATAAGAGGAGGTCAACTCGGTCTCCTTGTAGGGCAGTGGAACTCTCCTTTATACCCAGAGAGTAAATGGGCACCAGTATCACATAAGGGAGATGAGCTCTTCTATGGGTACAGTGTTTTCATGGGGTAGCTGTTTGTGGTCTCATCACAGAGCATTGGTCTTTTATGTTGTTTTGAGGATACCCAGGTAGCACAGTTGTTGAGTCTAGTGTCTAAAGGAGCAGAAACGTTATATGGTCTGTGAACCAGTATCATTCAGATGATGATGTCTTCTAGATGTTGAGCTCTGCATGTATGTAAACCCAGCCTTTGAGTAATGCCTCTTCCTTAGGCTAGTTTTGGTTTTACCCTAGAAGTTGGAAGGGTGCGTAGATATTTTTTTAACTGAAATTTGTTCTTAAATACAAAAAGATAACATAGAATTAAGCTGATCTTTAGGACAGATAACATTGATGCAAAACTGTGAAGTA
>abh_1627_B
CAnGAAGTAGACTAGCTACAGGCTCAATTTGATAAGAGGAGGTCAACTCGGTCTCCTTGTAGGGCAGTGGAACTCTCCTTTATACCCAGAGAGTAAATGGGCACCAGTATCACGTAAGGGAGATGAGCTCTTTTATGGGTACAGTGTATTCATGGGGTAGCTGTTTGTGGTCTCATCACAGAGCTTTGGTCTTTTATGTTGTTTTGAGGATACCCAGGTAGCACAGTTGTTGAGTCTAGTGTCTAAAGGAGCAGAAACGTTATATGGTCTGTGAACCAGTATCATTCAGATGATGATGTCTTCTAGATGTTGAGCTCTGCATGTATGTAAACCCAGCCTTTGAGTAATGCCTCTTCCTTAGGCTAGTTTTGGTTCTACCCTAGAAGTTGGAAGGGTGCGTAGATATTTTTTTAACTGAAATTTGTTCTTAAATACAAAAAGATAACATAGAATTAAGCTGATCTTTAGGACAGATAACATTGATGCAAAACTGTGAAGTA
```

- file with statistics (called *gene\_name\_STAT.stats*): a single file is created for each analyzed gene; each file contains the statistics for each analyzed individual. The counter for `lowercases` is for not phased sited with the ones with low quality.  
Example (file name: *abh\_STAT.stats*):  
```
************************************
Individual number: 1622
************************************

Length: 500
Heterozygotic sites: 2
Phased sites: 3
Not phased sites: 0
Two alternatives: 0
Low quality sites: 1
Lowercases: 1
Deletions: 0



************************************
Individual number: 1627
************************************

Length: 500
Heterozygotic sites: 2
Phased sites: 5
Not phased sites: 0
Two alternatives: 0
Low quality sites: 1
Lowercases: 1
Deletions: 0
```

- a file called *invalid\_files.txt*, containing all the file names in vcf, which were not correctly created or with sequence length less than 1 nucleotide
