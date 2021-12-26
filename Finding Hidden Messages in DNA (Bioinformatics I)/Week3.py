#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 20:06:55 2021

@author: yuzhu
"""

import math
positive_infnity = math.inf

# data=[]
# with open("dataset12.txt","r") as f:
#     data = f.readlines()
# for i in range(len(data)):
#     data[i] = data[i].strip('\n')

def remove(string):
    return " ".join(string.split())


def FindGenome(data):
    Genome = []
    for i in range(1, len(data)):
        Genome.append(remove(data[i]))
    return Genome


#Dna = FindGenome(data)
#print(Dna)
        
# file = open("dataset6.txt","r")
# for line in file:
#     fields = line.split(" ")
# #print(fields)

def HammingDistance(p,q):
    if len(p) != len(q):
        raise 'The length of strings are incorrect!'
    else:
        Mismatch = 0
        for i in range(len(p)):
            if p[i] != q[i]:
                Mismatch += 1
    return Mismatch

def Neighbors(Pattern, d):
    if d == 0:
        return [Pattern]
    elif len(Pattern) == 1:
        return ['A', 'C', 'G', 'T']
    else:
        Neighborhood = []
        SuffixNeighbors = Neighbors(Pattern[1::], d)
        for string in SuffixNeighbors:
            if HammingDistance(Pattern[1::], string) < int(d):
                Nucleotides = ['A','T','G','C']
                for x in Nucleotides:
                    Neighborhood.append(str(x)+string)
            else:
                Neighborhood.append(Pattern[0]+string)
        return Neighborhood

def MotifEnumeration(Genome, k, d):
        Patterns = []
        for j in range(len(Genome)):
            for i in range(len(Genome[j])-int(k)+1):
                kmer = Genome[j][i:i+k]
                neighbors = Neighbors(kmer,d)
                #print('neighbors of '+ kmer+ ' are:' + str(neighbors))
                for neighbor in neighbors:
                    GenomeCopy = Genome.copy()
                    Neineighbors = Neighbors(neighbor,d)
                    for i in range(len(Genome),0,-1):
                        for Neineighbor in Neineighbors:
                            if Neineighbor in Genome[i-1]:
                                del(GenomeCopy[i-1])
                                break
                    if GenomeCopy == []:
                        if neighbor not in Patterns:
                            Patterns.append(neighbor)
        return Patterns
    
#Genome = ['ACGT','ACGT','ACGT']
# pattern = MotifEnumeration(Genome, 5, 2)
# print(*pattern, sep =' ')
        
Motifs = [
"TCGGGGGTTTTT",
"CCGGTGACTTAC",
"ACGGGGATTTTC",
"TTGGGGACTTTT",
"AAGGGGACTTCC",
"TTGGGGACTTCC",
"TCGGGGATTCAT",
"TCGGGGATTCCT",
"TAGGGGAACTAC",
"TCGGGTATAACC"
]

def CountEntropy(Motifs):
    Column = []
    for i in range(len(Motifs[0])):
        List = [0,0,0,0]
        n = len(Motifs)
        for j in range(len(Motifs)):
            if Motifs[j][i] == 'A':
                List[0] += (1 / n)
            elif Motifs[j][i] == 'T':
                List[1] += (1 / n)
            elif Motifs[j][i] == 'C':
                List[2] += (1 / n)
            else:
                List[3] +=  (1 / n)  
        Column.append(List)
    Score = 0
    for i in range(len(Column)):
        score = 0
        column = Column[i]
        for j in range(4):
            if column[j] != 0:
                num = -column[j] * math.log2(column[j])
                score += num
        Score += score
    return Score

#print(CountEntropy(Motifs))
 
def DistanceBetweenPatternAndStrings(Pattern, Dna):
    k = len(Pattern)
    distance = 0
    for string in Dna:
        for i in range(len(string)-k+1):
            Text = string[i:i+k]
            if i == 0:
                Hammingdistance = HammingDistance(Pattern, Text)
            if Hammingdistance > HammingDistance(Pattern, Text):
                Hammingdistance = HammingDistance(Pattern, Text)
        distance = distance + Hammingdistance
    return distance    

#Dna = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
#Pattern = 'AAA'
#print(DistanceBetweenPatternAndStrings('AGAAG', fields))
       

def ConvertDecimalToBase(Decimal, Base):
    Result = []
    while Decimal > 0:
        Number = Decimal % Base
        Result.append(Number)
        Decimal = (Decimal - Number) // Base
    str1 = ''
    for i in range(len(Result)):
        str1 += str(Result[-i-1])
    return str1

#print(ConvertDecimalToBase(25, 4))
    
def AllStrings(k):
    Strings = []
    for i in range(4**k):
        combo = ''
        n = str(ConvertDecimalToBase(i, 4))
        if len(n) != k:
            str1 = ''
            for m in range(k-len(n)):
                str1 += '0'
            n = str1 + n
        for j in range(1,k+1):
            if n[-j] == '0':
                combo += 'A'
            elif n[-j] == '1':
                combo += 'T'
            elif n[-j] == '2':
                combo += 'G'
            else:
                combo += 'C'
        Strings.append(combo)
    return Strings
                
#print(AllStrings(3))
#print(len(AllStrings(3)))

#Dna = ['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG','GCTGAGCACCGG','AGTTCGGGACAG']
#Dna = ['ACGT','ACGT','ACGT']
#Dna = ['ATA','ACA','AGA','AAT','AAC']
#Dna = ['AAG','AAT']

def MedianString(Dna, k):
    Patterns = AllStrings(k)
    for i in range(len(Patterns)):
        Pattern = Patterns[i]
        if i == 0:
            distance = DistanceBetweenPatternAndStrings(Pattern, Dna)
            Median = Pattern
        elif distance > DistanceBetweenPatternAndStrings(Pattern, Dna):
            distance = DistanceBetweenPatternAndStrings(Pattern, Dna)
            Median = Pattern
    return Median

def AllMedianString(Dna, k):
    Medians = {}
    Patterns = AllStrings(k)
    for i in range(len(Patterns)):
        Pattern = Patterns[i]
        Medians[Pattern] = DistanceBetweenPatternAndStrings(Pattern, Dna)
    values = Medians.values()
    min1 = min(values)
    Median = []
    for key in Medians:
        if Medians[key] == min1:
            Median.append(key)
    return Median


Dna = ['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG']
print(AllMedianString(Dna, 7))
print(Dna)
    
profile = {
    'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    'T': [0.1, 0.2, 0.1, 0.1, 0.2]
}

def ProfileMostProbableKmer(Text, k, profile):
    for i in range(len(Text)-int(k)+1):
        Pattern = Text[i:i+int(k)]
        Number = 1
        for j in range(int(k)):
            base = Pattern[j]
            Number *= profile[base][j]
        if i == 0:
            Max = Number
            pattern = Pattern
        elif Max < Number:
            Max = Number
            pattern = Pattern
        #print(pattern, Pattern, Number, Max)
    return pattern

def convertToList(string):
    List = string.split(' ')
    Copy = []
    for item in List:
        Copy.append(float(item))
    return Copy

#print(ProfileMostProbableKmer('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, profile))
# print(data)
# Text = data[0]
# k = data[1]
# profile = {}
# profile['A'] = convertToList(data[2])
# profile['C'] = convertToList(data[3])
# profile['G'] = convertToList(data[4])
# profile['T'] = convertToList(data[5])
# print(Text,k, profile)
# print(ProfileMostProbableKmer(Text, k, profile))
  
def formProfile(Dna):
    profile = {}
    profile['A'] = []
    profile['C'] = []
    profile['G'] = []
    profile['T'] = []
    for i in range(len(Dna[0])):
        (A,C,G,T) = (1,1,1,1)
        for j in range(len(Dna)):
            if Dna[j][i] == 'A':
                A += 1
            elif Dna[j][i] == 'C':
                C += 1
            elif Dna[j][i] == 'G':
                G += 1
            else:
                T += 1
        Total = A+C+G+T
        profile['A'].append(A/Total)
        profile['C'].append(C/Total)
        profile['G'].append(G/Total)
        profile['T'].append(T/Total)
    return profile

#print(formProfile(Dna))

def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(len(Dna)):
        BestMotifs.append(Dna[i][0:k])
    for i in range(len(Dna[0])-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1,t):
            profile = formProfile(Motifs)
            Motif = ProfileMostProbableKmer(Dna[j], k, profile)
            Motifs.append(Motif)
        if CountEntropy(Motifs) < CountEntropy(BestMotifs):
            print(CountEntropy(Motifs),Motifs[0])
            BestMotifs = Motifs
    return BestMotifs

# Dna1 = ['GGCGTTCAGGCA',
# 'AAGAATCAGTCA',
# 'CAAGGAGTTCGC',
# 'CACGTCAATCAC',
# 'CAATAATATTCG']
# Dna2 = ['GCCCAA',
# 'GGCCTG',
# 'AACCTA',
# 'TTCCTT']  
# Dna3 = ['GCAGGTTAATACCGCGGATCAGCTGAGAAACCGGAATGTGCGT',
# 'CCTGCATGCCCGGTTTGAGGAACATCAGCGAAGAACTGTGCGT',
# 'GCGCCAGTAACCCGTGCCAGTCAGGTTAATGGCAGTAACATTT',
# 'AACCCGTGCCAGTCAGGTTAATGGCAGTAACATTTATGCCTTC',
# 'ATGCCTTCCGCGCCAATTGTTCGTATCGTCGCCACTTCGAGTG']
# Dna4 = ['GAGGCGCACATCATTATCGATAACGATTCGCCGCATTGCC',
# 'TCATCGAATCCGATAACTGACACCTGCTCTGGCACCGCTC',
# 'TCGGCGGTATAGCCAGAAAGCGTAGTGCCAATAATTTCCT',
# 'GAGTCGTGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG',
# 'GACGGCAACTACGGTTACAACGCAGCAACCGAAGAATATT',
# 'TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT',
# 'AAGCGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG',
# 'AATTGAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA']
# Dna = ['GACCTACGGTTACAACGCAGCAACCGAAGAATATTGGCAA',
# 'TCATTATCGATAACGATTCGCCGGAGGCCATTGCCGCACA',
# 'GGAGTCTGGTGAAGTGTGGGTTATGGGGCAGACTGGGAAA',
# 'GAATCCGATAACTGACACCTGCTCTGGCACCGCTCTCATC',
# 'AAGCGCGTAGGCGCGGCTTGGCATCTCGGTGTGTGGCCAA',
# 'AATTGAAAGGCGCATCTTACTCTTTTCGCTTAAAATCAAA',
# 'GGTATAGCCAGAAAGCGTAGTTAATTTCGGCTCCTGCCAA',
# 'TCTGTTGTTGCTAACACCGTTAAAGGCGGCGACGGCAACT']

# Dna = ['GGCGTTCAGGCA',
# 'AAGAATCAGTCA',
# 'CAAGGAGTTCGC',
# 'CACGTCAATCAC',
# 'CAATAATATTCG']

# Dna = data[1::]
# print(Dna)
# result = GreedyMotifSearch(Dna, 12, 25)
# print(*result, sep=' ')
