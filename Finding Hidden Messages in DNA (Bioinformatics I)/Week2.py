#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 14:22:20 2021

@author: yuzhu
"""


data=[]
with open("Salmonella_enterica.txt","r") as f:
    data = f.readlines()
    
def remove(string):
    return "".join(string.split())

#print(data[3])

def FindGenome(data):
    Genome = ''
    i = len(data)
    for i in range(1,i):
        Genome += str(data[i])
    return "".join(Genome.split())
#p = remove(data[0])
#q = remove(data[1])
# print(p)
# print(q)
#print(FindGenome(data))
Genome = FindGenome(data)

i = 'AFEDS'


def SkewValue(Genome):
    Skew = 0
    SkewValues = [0]
    for i in range(len(Genome)):
        if Genome[i] == 'C':
            Skew -= 1
        elif Genome[i] == 'G':
            Skew += 1
        SkewValues.append(Skew)
    return SkewValues

# Skew = SkewValue('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
# print(*Skew, sep=' ')
# print(Skew[11])
# print(Skew[24])

def MinimumSkew(Genome):
    Skew = SkewValue(Genome)
    Min = min(Skew)
    Position = []
    for i in range(len(Skew)):
        if Skew[i] == Min:
            Position.append(i)
    return Position 
 
#Genome = 'GCATACACTTCCCAGTAGGTACTG'
#Genome = 'TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'
Position = MinimumSkew(Genome)
print(Position)

      
def HammingDistance(p,q):
    if len(p) != len(q):
        raise 'The length of strings are incorrect!'
    else:
        Mismatch = 0
        for i in range(len(p)):
            if p[i] != q[i]:
                Mismatch += 1
    return Mismatch

#p = 'CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG'
#q = 'ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT'    
#p = 'GGGCCGTTGGT'
#q = 'GGACCGTTGAC'
#print(HammingDistance(p,q))
    
def ApproximatePatternMatching(Pattern,Genome,d):
    m = len(Pattern)
    n = len(Genome)
    Position = []
    for i in range(n-m+1):
        if HammingDistance(Genome[i:i+m],Pattern) <= int(d):
            Position.append(i)
    return Position

# Pattern = remove(data[0])
# Genome = remove(data[1])
# d = remove(data[2])
# num = ApproximatePatternMatching(Pattern,Genome,d)
# print(*num, sep=' ')

def FrequentWords(Genome,Pattern,d):
    k = len(Pattern)
    n = len(Genome)
    Words = []
    for i in range(n-k+1):
        g = Genome[i:i+k]
        if HammingDistance(Pattern,g) <= int(d):
            Words.append(g)
    return Words

#print(FrequentWords('CATGCCATTCGCATTGTCCCAGTGA', 'CCC',2))
 
    
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

# neigh = Neighbors('ACTACCAG',3)
# print(*neigh,sep=' ')
        
def ReverseComplement(Pattern):
    PatternRv = ''
    for i in range(len(Pattern)):
        if Pattern[i] == 'A':
            PatternRv += 'T'
        elif Pattern[i] == 'T':
            PatternRv += 'A'
        elif Pattern[i] == 'C':
            PatternRv += 'G'
        else:
            PatternRv += 'C'
    return PatternRv[::-1]

def FrequentWordsWithMismatches(Genome, k, d):
    Patterns = []
    freqMap = {}
    n = len(Genome)
    for i in range(n-k+1):
        Pattern = Genome[i:i+k]
        neighborhood = Neighbors(Pattern, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] = freqMap[neighbor] + 1
    m = max(freqMap.values())
    for Pattern in freqMap:
        if freqMap[Pattern] == m:
            Patterns.append(Pattern)
    return Patterns 
    
# freq = FrequentWordsWithMismatches('ACGT', 4, 3)
# print(*freq, sep=' ')    
# print(len(freq))      
    
def FrequentWordsWithMismatchesAndComplement(Genome,k,d):
    Patterns = []
    freqMap = {}
    n = len(Genome)
    for i in range(n-k+1):
        Pattern = Genome[i:i+k]
        neighborhood = Neighbors(Pattern, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] = freqMap[neighbor] + 1
        PatternRc = ReverseComplement(Pattern)
        neighborhood2 = Neighbors(PatternRc,d)
        for m in range(len(neighborhood2)):
            neighbor = neighborhood2[m]
            if neighbor not in freqMap:
                freqMap[neighbor] = 1
            else:
                freqMap[neighbor] = freqMap[neighbor] + 1
    m = max(freqMap.values())
    for Pattern in freqMap:
        if freqMap[Pattern] == m:
            Patterns.append(Pattern)
    return Patterns 
# print(Genome)
# freq = FrequentWordsWithMismatchesAndComplement(Genome,5,2)
# print(*freq, sep=' ')

def DefineOrigin(Position) :
    dic = {}
    for item in Position:
        Sequence = Genome[3764856-500:3764856+500]
        for i in range (5,10):
            for j in range(6):
                freq = FrequentWordsWithMismatchesAndComplement(Sequence,i,j)
                for k in range(len(freq)):
                    if freq[k] not in dic:
                        dic[freq[k]] = 1
                    else:
                        dic[freq[k]] += 1
    Max = max(dic.values())
    Patterns = []
    for item in dic:
        if dic[item] == Max:
            Patterns.append(item)
    return Patterns
        
seq = DefineOrigin(Position)
print(seq)
    