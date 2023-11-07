#!/usr/bin/env python
import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Take uniquely mapped SAM file, remove PCR duplicates and output: (1) SAM file with only unique reads, (1) duplicate.sam file containing other PCR duplicates, and (1) report.tsv file ")
    parser.add_argument("-f", "--file", help="Path to input SAM file", type= str, required=True)
    parser.add_argument("-u", "--umi", help="Path to umi file", type=str)
    parser.add_argument("-o", "--outfile", help="Path to output directory", type=str)
    return parser.parse_args()
args=get_args()	

file_name = args.file  
umi_file = args.umi
output_file = args.outfile

def evaluateStrand (flag) -> str():
    ''' This function takes in the bitwise flag and inform what strand this is'''
    strand = str()
    if ((flag & 16) != 16):
        strand = "forward"
    else:
        strand = "reverse"
    return strand

# Create a function that will adjust starting position take into account soft clipping and strandedness
def adjust_pos (left_most_pos, cigar, flag) -> int():
    '''This function takes in left most position (extracted from SAM file) + cigar string + bitwise flag
    and caclulate the true starting position of the read '''
    add = 0
    soft_clip = 0
    pattern = re.compile(r'[0-9]+[MINDSHP=X]')
    # Each component of cigar string is stored in a list (ex: ["1S","2M","6D","2S"])
    cigar_list = pattern.findall(cigar)

    # evaluate strandedness
    if ((flag & 16) != 16):
        revFlag = False
    else:
        revFlag = True
    
    # if reverse string then flip the list around so if there is soft clipping we'll use the right-hand side (3') S-value
    if revFlag == True:
        cigar_list = cigar_list [::-1]

    for item in cigar_list:
        letter = item[-1]
        number = int(item[:-1])
        if letter == "N" or letter == "M" or letter == "D":
            add += number
    if "S" in cigar_list[0]:
        soft_clip = int(cigar_list[0][:-1])

    if revFlag == False:
        adjusted_start = left_most_pos - soft_clip
    elif revFlag == True:
        adjusted_start = left_most_pos + add + soft_clip
    return adjusted_start

# Store all the known UMIs into a set
known_UMI = set()
with open(umi_file,"r") as file:
    while True:
        line = file.readline().strip()
        if line == "":
            break 
        known_UMI.add(line)

# Initiate some variables:
non_duplicate = set()
current_chr = str()

# Initiate some counters:
records = int()
unknownUMI_counter = int()
duplicate_counter = int()
unique_counter = int()
chromosome_counter = {}

# Main block of code:
# This script will output a SAM file with the unique reads and another SAM file with other PCR duplicates. 
with open (args.file,"r") as f_in, open(f'{output_file}/uniqueReads.sam',"w") as f_out, open (f'{output_file}/Duplicates.sam',"w") as f_dupe:
    for line in f_in:
        line = line.strip("\n")
        if line.startswith("@"):
            f_out.write(f'{line}\n')
            continue
        else:
        # These are all the components in the SAM file
            components = line.split("\t")
            UMI = components[0].split(":")[-1]
            flag = int(components[1])
            chromosome = components[2]
            left_most_pos = int(components[3])
            cigar = components[5]

        # Use the evaluateStrand function to evaluate strandedness ("forward" or "reverse")
        # Use adjusted_startPos function to obtain the true starting position
            adjusted_startPos = adjust_pos(left_most_pos, cigar, flag)
            strand = evaluateStrand (flag)
        # Store all the identifiers: UMI, adjusted starting position, chromosome, and strand in a tuple for each iteration
            read = (UMI,adjusted_startPos,chromosome,strand)

        # This chunk of code reset the set of all unique read identity whenever we move on to a new chromosome
            if chromosome != current_chr:
                current_chr = chromosome
                non_duplicate.clear()
        
        # Write unique reads into a file:
        # Unique reads must has known UMIs with identfiers not already existed in the non_duplicate set
            if UMI in known_UMI:
                if read in non_duplicate:
                    duplicate_counter += 1
                    f_dupe.write(f'{line}\n')
                    records += 1
        # Write duplicates into a file:
                else:
                    f_out.write(f'{line}\n')
                    non_duplicate.add(read)
                    unique_counter += 1
                    records += 1
                    if not chromosome_counter.get(chromosome):
                        chromosome_counter[chromosome] = 1
                    else:
                        chromosome_counter[chromosome] += 1
            else:
                unknownUMI_counter += 1
                records += 1

# Write a summary report file:          
with open (f'{output_file}/report.tsv',"w") as report:
    report.write("Statistic Summary\n")
    report.write(f'Total number of reads: {records}\n')
    report.write(f'Number of unique reads: {unique_counter}\n')
    report.write(f'Number of PCR duplicates: {duplicate_counter}\n')
    report.write(f'Unknown UMIs: {unknownUMI_counter}\n')

    report.write("Chromosome\tReads per chromosome\n")
    for chr,count in chromosome_counter.items():
        report.write(f'{chr}\t{count}\n')