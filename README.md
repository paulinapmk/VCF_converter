# VCF Converter to FASTA format with some statistics
Next Generation Sequencing (NGS) generates a huge amount of data. One of the formats that allows the storage of information about sequences and their annotations is Variant Call Format (VCF). This thesis describes the main features of VCF file format, in particular focusing on the description of information about polymorphic sites in sequences. They have an important role in population studies, among other things. The key element of this dissertation is a converter from VCF to FASTA format. It is a script written in Python, which can be executed from command line. Therefore it is possible to execute the program in different pipelines. The appended manual includes the description of all available options as well as the explanation of the converter’s functionalities.

## General information
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

For the correct execution of the script, the header line starting with `\#CHROM` has to be present. `1549` in the example is the id of my individual (its presence is important too).

By default, if the difference in DP value between two neighbour nucleotides is greater than 80%, the letter `n` is transmitted to the output file. This value can be modified in the script according to your required precision.
