__author__ = 'Paulina Kania'

import os
import parser

options = parser.parser()
input_type = options.input_type
folder_path = options.input_path
file_path = options.filename
output_path_def = options.output_path

if folder_path == '':
    path = os.listdir(os.getcwd())
else:
    path = os.listdir(folder_path)

if output_path_def == '':
    output_path = os.listdir(os.getcwd())
else:
    output_path = os.listdir(output_path_def)

splited = []
individual_ids = []
invalids = []
def input_split(input):
    """
    This function splits the file skipping all the meta-information lines
    :param input: file in vcf format
    :return: file lines splitted into arrays of strings
    """
    for line in input:
        if "#CHROM" in line:
            splited.append(line.split())
        elif '#' in line:
            continue
        else:
            splited.append(line.split())
    return splited

def deletion_search(splited_input):
    """
    This function returns the missing positions in the input file
    :param splited_input: splited input file
    :return: array with missing indexes
    """
    a = []
    for line in splited_input[1:]:
        a.append(int(line[1]))
    print [(e1+1) for e1,e2 in zip(a, a[1:]) if e2-e1 != 1]

def getPhases(splited_input, index):
    """
    This is the function with the main algorithm for reading the sequence
    :param splited_input: file lines splitted into arrays of strings
    :param index: column to be read (depends on the individual position in file)
    :return: the two sequences for each individual and the statistics
    """
    outA, outB, stats = "",  "", ""
    counter = 0
    real_pos = 0
    heterozygotic_sites, phased_sites, not_phased_sites = 0, 0, 0
    two_alt_counter, low_quals_counter, lowercase_counter = 0, 0, 0
    for line in splited_input[1:]:
        pos = int(float(line[1]))
        real_pos += 1
        ref = line[3]
        alt = line[4].split(',')
        filter = line[6]
        genotypes = line[9:]
        gt = genotypes[index].split(':')[0]
        format = line[8]
        dp = 0
        if 'GQ' in format:
            gq_index = format.split(':').index('GQ')
            gq = int(genotypes[index].split(':')[gq_index])
        if 'DP' in format:
            dp_index = format.split(':').index('DP')
            try:
                dp = int(genotypes[index].split(':')[dp_index])
            except IndexError:
                dp = 0
            if line == splited_input[1]:
                last_dp = dp
        while real_pos < pos:#checking if deletions are present
            outA += "-"
            outB += "-"
            real_pos +=1
            last_dp = 0
        if ('GQ' in format) and (gq <= 40):
            outA += "n"
            outB += "n"
            lowercase_counter += 1
            low_quals_counter += 1
        else:
            if filter == 'LowQual':
                outA += "n"
                outB += "n"
                lowercase_counter += 1
                low_quals_counter += 1
            else:
                if (('DP' in format) and (dp < (last_dp - last_dp*80/100))) or \
                        (('DP' in format) and (dp < 5)):#if dp decreased by 80% or less than 5 (5 is the minimal acceptable DP)
                    outA += "n"
                    outB += "n"
                    lowercase_counter += 1
                    low_quals_counter += 1

                else:
                    last_dp = dp
                    if gt == '0|0':
                        outA += ref
                        outB += ref
                    elif gt == '1|1':
                        outA += alt[0]
                        outB += alt[0]
                    elif gt == '0|1':
                        outA += ref
                        outB += alt[0]
                        phased_sites += 1
                    elif gt == '1|0':
                        outA += alt[0]
                        outB += ref
                        heterozygotic_sites += 1
                        phased_sites += 1
                    elif gt == '2|2':
                        outA += alt[1]
                        outB += alt[1]
                    elif gt == '0/2':
                        outA += ref.lower()
                        outB += alt[1].lower()
                        lowercase_counter += 1
                        two_alt_counter += 1
                        heterozygotic_sites += 1
                        not_phased_sites += 1
                    elif gt == '1/2':
                        outA += alt[0].lower()
                        outB += alt[1].lower()
                        lowercase_counter += 1
                        two_alt_counter += 1
                        heterozygotic_sites += 1
                        not_phased_sites += 1
                    elif gt == '2/2':
                        outA += alt[1].lower()
                        outB += alt[1].lower()
                        lowercase_counter += 1
                    elif gt == '1/1':
                        outA += alt[0].lower()
                        outB += alt[0].lower()
                        lowercase_counter += 1
                    elif gt == '0/0':
                        outA += ref.lower()
                        outB += ref.lower()
                        lowercase_counter += 1
                    elif gt == '0/1':
                        if counter == 0:
                            outA += ref
                            outB += alt[0]
                            counter += 1
                            heterozygotic_sites += 1
                            phased_sites += 1
                        else: # if not phased for other reason than because it is the first
                            outA += ref.lower()
                            outB += alt[0].lower()
                            lowercase_counter += 1
                            not_phased_sites += 1
                    elif gt == './.':
                        outA += '-'
                        outB += '-'
                    else:
                        print 'Unrecognized genotype info at position ', pos
    deletions = outA.count('-', 0, len(outA))
    stats_out = [len(outA), heterozygotic_sites, phased_sites, not_phased_sites, two_alt_counter, low_quals_counter, lowercase_counter, deletions]
    stats = stats_output(stats_out)
    return outA, outB, stats

def stats_output(stats):
    """
    This function just prepares the statistics' appearance to be written into the output file
    :param stats: array of all the elements to be put into the statistics
    :return: string containing the statistics with description
    """
    stats = "Length: " + str(stats[0]) + "\n" + \
            "Heterozygotic sites: " + str(stats[1]) + "\n" + \
            "Phased sites: " + str(stats[2]) + "\n" + \
            "Not phased sites: " + str(stats[3]) + "\n" + \
            "Two alternatives: " + str(stats[4]) + "\n" + \
            "Low quality sites: " + str(stats[5]) + "\n" + \
            "Lowercases: " + str(stats[6]) + "\n" + \
            "Deletions: " + str(stats[7]) + "\n\n\n"
    return stats

def output_definition(splited_input):
    """
    This function runs the main algorithm, checks if in the chosen output directory the output files already exist
    for the analyzed chromosome and individual. If not it creates it
    :param splited_input: file lines splitted into arrays of strings
    :return: it just returns the output files written in the previously chosen output directory
    """
    individual_ids = splited_input[0][9:]
    outA, outB, name_A, name_B, output, stats_out, output_stats, name_stats = "", "", "", "", "", "", "", ""
    for ind_num in individual_ids[0:]:
        ind_index = individual_ids.index(ind_num)
        outA, outB, stats_out = getPhases(splited_input, ind_index)
        for line in splited_input[1:]:
            gene = line[0]
            output = output_path_def + str(gene) + ".fasta"
            output_stats = output_path_def + str(gene) + "_STAT.stats"
            name_A = ">"+gene+"_"+str(individual_ids[individual_ids.index(ind_num)]+"_A")
            name_B = ">"+gene+"_"+str(individual_ids[individual_ids.index(ind_num)]+"_B")
            name_stats = "************************************\n" + \
                         "Individual number: " + str(individual_ids[individual_ids.index(ind_num)]) + "\n" + \
                         "************************************\n"
        if os.path.exists(output):
            with open(output,'a+') as f:
                if any(name_A == x.rstrip("\r\n") for x in f) and any(name_B == x.rstrip("\r\n") for x in f):
                    print "Already existing outputs: " + name_A[1:] + ' , ' + name_B[1:]
                else:
                    f.write(name_A + "\n")
                    f.write(outA + "\n")
                    f.write(name_B + "\n")
                    f.write(outB + "\n")
        else:
            with open(output,'w') as f:
                f.write(name_A + "\n")
                f.write(outA + "\n")
                f.write(name_B + "\n")
                f.write(outB + "\n")

        if os.path.exists(output_stats):
            with open(output_stats,'a+') as f:
                if any(("Individual number: " + ind_num) == x.rstrip("\r\n") for x in f):
                    print "Already existing statistics' output for individual: " + ind_num
                else:
                    f.write(name_stats + "\n")
                    f.write(stats_out + "\n")
        else:
            with open(output_stats,'w') as f:
                f.write(name_stats + "\n")
                f.write(stats_out + "\n")
    del splited[:] #clear the splited array

def input_definition():
    """
    It runs all the other functions basing on the choices made by the user
    """
    if input_type == 'file':
        if file_path.endswith(".vcf"):
            input = open(os.path.join(folder_path, file_path))
            split = input_split(input)
            if len(split) <= 1:
                print "Invalid file: ", file_path
                return
            output_definition(split)
            input.close()

        return
    else:
        for filename in path:
            if filename.endswith(".vcf"):
                #print "Reading " + filename
                input = open(os.path.join(folder_path, filename))
                split = input_split(input)
                if len(split) <= 1:
                    print "Invalid file: ", filename
                    invalids.append(filename)
                    del splited[:] #clear the splited array
                    continue
                else:
                    output_definition(split)
                    input.close()
            else:
                continue
        if len(invalids) != 0:
            with open("invalid_files.txt",'w') as f:
                f.write("Invalid files: \n\n")
                f.writelines( "%s\n" % item for item in invalids )
        return

if __name__ == "__main__":
    input_definition()