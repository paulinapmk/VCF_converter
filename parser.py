__author__ = 'Paulina Kania'

from optparse import OptionParser

def parser():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.description = "This is a converter from VCF format to FASTA."
    parser.add_option("-t", "--type", dest="input_type",
                      help="chose 'folder' or 'file' to be executed")
    parser.add_option("-f", "--file", dest="filename",
                      help="read data from FILENAME", default="")
    parser.add_option("-i", "--input",
                      dest="input_path", help="read input path", default="")
    parser.add_option("-o", "--output",
                      dest="output_path", help="read output path", default="")
    (options, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("wrong number of arguments")
    if options.input_type == None:
        parser.print_help()
        parser.error("you have to specify input type")
    return options

if __name__ == "__main__":
    parser()
