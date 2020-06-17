#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This algorithm is capable of processing a fasta file and generating
Illumina-like reads from reference sequences. Additionally, it allows
modifying the percentage of identity of the generated sequences.'''

__author__ = "Ary Rivillas"
__copyright__ = "Copyright 2020, Biotecnologia Microbiana Research Group,\
     Universidad Nacional de Colombia - Sede Medellín"
__credits__ = ["Pablo Gutierrez", "Daniel Tejada", "Andrea Restrepo",\
     "Susana Giraldo"]
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Ary Rivillas"
__email__ = "amrivillast@unal.edu.co"
__status__ = "Development"
__date__ = '2020/06/5'

import re
import os
import random
from random import randint
from random import choices

# The options of the following block can be modified depending on what
# is required

##############################################################################
############################### MODIFY ME ####################################
##############################################################################

# Name of the file from which the reads will be generated
file_base = "viruses_of_resolution.fasta"

# Quantity of reads for every sequence in fasta file.
reads = 20000

# FASTQ files are generate by identity. The amount of files is equal to
# end_identity - start_identity
start_identity = 95
end_identity = 96

##############################################################################
##############################################################################
##############################################################################

dir_path = os.path.dirname(os.path.realpath(__file__))
path_to_file = os.path.join(dir_path, file_base)
count = len(open(path_to_file).readlines())

title = []
seq = []
seq_string = ""
z = 0
w = 0

# The file is open. Sequence and title are save in lists.
with open(path_to_file) as lines:

    for i, line in enumerate(lines):

        seq_re = '[^>]'
        seq_hit = re.match(seq_re, line)
      
        if seq_hit:
     
            line_fix2 = line.strip()
            seq_string += str(line_fix2)
         
            if i == (count - 1):

                seq.append(seq_string)

        else:

            line_fix = line.strip()
            title.append(line_fix)
         
            if (len(title) - len(seq)) == 2:

                seq.append(seq_string)
                seq_string = ""

# Possible fragment lengths and their associated probability, simulating
# a normal distribution
fragment = [225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500]
prob_fragment = [0.01, 0.02, 0.03, 0.06, 0.12, 0.26, 0.26, 0.12, 0.06,\
     0.03, 0.02, 0.01]

# Possible quality values and their associated probability, giving more
# weight to high values
quality = [33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
     49, 50, 51, 52,53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
     67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
     85, 86, 87, 88, 89, 90]
prob_quality = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
     0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.015,
     0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.0238, 0.0238,
     0.0238, 0.0238, 0.0238, 0.0238, 0.0238, 0.05, 0.05, 0.05, 0.05, 0.05,
     0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.0238, 0.0238, 0.0238, 0.0238,
     0.0238, 0.0238, 0.0238, 0.0238, 0.0238]

# This cycle determines how many files are generated
for ki in range(start_identity, end_identity, 1):

    for i in range(len(seq)):

        reads_with_reverse = reads*2
        sequence = seq[i]
        fq_s =""
        rq_s = ""
        answer = []
        z = 0
        xy = 1        
        complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        percent = ki    
        percent_fix = 100 - int(percent)        
        length = len(sequence)
        perc = round((int(percent_fix)*length)/100)
        
        for i in range(perc + 1):    
    
            seq_rand = randint(0, length - 1)
            seq_nt_rand = sequence[seq_rand:seq_rand + 1]            
            start_ = sequence[0:seq_rand]            
            end_ = sequence[seq_rand + 1:len(sequence) + 1]            
            nucleotides = ["A", "C", "G", "T"]
            nt_rand = random.randint(0, 3)
            lucky_nt = nucleotides[nt_rand]            
            seq_replace = seq_nt_rand.replace(seq_nt_rand, lucky_nt)
            
            if seq_replace == seq_nt_rand:
                
                nucleotides_fix = nucleotides
                nucleotides_fix.remove(seq_nt_rand)               
                nt_rand_if = random.randint(0, 2)                
                lucky_nt_if = nucleotides_fix[nt_rand_if]
                seq_replace_if = seq_nt_rand.replace(seq_nt_rand, lucky_nt_if)                
                complete1 = start_ + seq_replace_if + end_
                sequence = complete1  
                
            else:

                complete2 = start_ + seq_replace + end_
                sequence = complete2
        
        while z < reads_with_reverse:

            inicio = randint(0, len(sequence))
            w = random.choices(fragment, prob_fragment)
            substr = sequence[inicio:(inicio + int(w[0]))]
                        
            for nucleotide in substr:

                answer.append(complement[nucleotide])
        
            substr_comp = ''.join(answer)
            read_f = substr[0:101]            
            reverso = substr_comp[::-1]
            read_r = reverso[0:101]        
        
            if len(read_f) == 101 and len(read_r) == 101:

                f_q = random.choices(quality, prob_quality, k = 101)               
                flowcell = randint(0, 9999)
                xcor = randint(0, 99999)
                ycor = randint(0, 99999)
               
                for i in f_q:

                    phred_f = chr(i)
                    fq_s += str(phred_f)

                file_base_F = "reads_F_" + str(ki) + ".fastq"
                path_to_file_F = os.path.join(dir_path, file_base_F)
               
                with open(path_to_file_F, "a+") as fastq_file:    
        
                    print("@MG00HS10:589:C69NNACXX:3:"
                    + str(flowcell)
                    + ":"
                    + str(xcor)
                    + ":"
                    + str(ycor)
                    + " :1:N:0:CTTGTA"
                    + '\n' + read_f + '\n' + "+" + '\n' 
                    + str(fq_s), file = fastq_file)
               
                r_q = random.choices(quality, prob_quality, k = 101)

                for j in r_q:

                    phred_r = chr(j)
                    rq_s += str(phred_r)

                file_base_R = "reads_R_" + str(ki) + ".fastq"
                path_to_file_R = os.path.join(dir_path, file_base_R)

                with open(path_to_file_R, "a+") as fastq_file2: 
                   
                    print("@MG00HS10:589:C69NNACXX:3:"
                    + str(flowcell)
                    + ":"
                    + str(xcor)
                    + ":"
                    + str(ycor)
                    + " :2:N:0:CTTGTA"
                    + '\n' + read_r + '\n' + "+" + '\n'
                    + str(fq_s), file = fastq_file2)

                z += 2
                xy += 1
                fq_s = ""
                rq_s = ""
                answer = []
