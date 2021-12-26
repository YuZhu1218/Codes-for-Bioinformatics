#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 03:37:37 2021

@author: yuzhu
"""

data=[]
with open("E. coli.txt","r") as f:
    data = f.readlines()
    
def remove(string):
    return "".join(string.split())

Genome=remove(data[0])
#Genome=remove(data[1])
#print(Text)
# print(Pattern)
#print(Genome)
# k=len(Pattern)
# print(k)

def PatternCount(Text, Pattern):
    k=len(Pattern)
    count= 0
    for i in range(len(Text)-k):
        pattern=Text[i:i+k]
        if Pattern==pattern:
            count+=1
    return count

# print(PatternCount('GACCATCAAAACTGATAAACTACTTAAAAATCAGT', 'AAA'))


def FrequentWords(Text, k):
    FrequentPatterns=[]
    maxCount = 0
    m = len(Text)-k
    for i in range(m):
        Pattern = Text[i:k+i]
        Count = PatternCount(Text, Pattern)
        # print('Pattern ' + Pattern + ' has a count of ' + str(Count)+ '; maxCount =' + str(maxCount))
        if Count > maxCount:
            FrequentPatterns = [Pattern]
            maxCount = Count
        elif Count == maxCount:
            if Pattern not in FrequentPatterns:
                FrequentPatterns.append(Pattern)
    return FrequentPatterns



def FrequencyTable(Text, k):
    freqMap={}
    n=len(Text)
    for i in range(n-k):
        Pattern = Text[i:k+i]
        if Pattern not in freqMap:
            freqMap[Pattern] = 1
        else:
           freqMap[Pattern] = freqMap[Pattern]+1
    return freqMap


def BetterFrequentWords(Text, k):
    FrequentPatterns = []
    freqMap = FrequencyTable(Text, k)
    allValues = freqMap.values()
    maxValue = max(allValues)
    # print(allValues)
    # print(maxValue)
    for Pattern in freqMap:
        if freqMap[Pattern] == maxValue:
           FrequentPatterns.append(Pattern) 
    return FrequentPatterns

# Text1 = 'CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT'
# print(BetterFrequentWords(Text1,3))

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

#print(ReverseComplement(Text))

def PatternMatching(Genome,Pattern):
    Position = []
    m = len(Genome)
    k = len(Pattern)
    for i in range(m-k):
        if Genome[i:i+k] == Pattern:
            Position.append(i)
    return Position

# data = PatternMatching(Genome,'CTTGATCAT')
# print(*data, sep='  ')

def ClumpFinding(Genome,k,L,t):
    Patterns = []
    n = len(Genome)
    for i in range(n-L):
        Window = Genome[i:i+L]
        freqMap = FrequencyTable(Window, k)
        for key in freqMap:
            if freqMap[key] >= t:
                if key not in Patterns:
                    Patterns.append(key)
    return len(Patterns)
        
#print(ClumpFinding(Genome,9,500,3))
    
def FindClumps(Genome,k,L,t):
    Patterns = []
    kmer = {}
    n = len(Genome)
    for i in range(n-k):
        pattern = Genome[i:i+k]
        if pattern not in kmer:
            kmer[pattern] = [i]
        else:
            kmer[pattern].append(i)
            if len(kmer[pattern]) == t:
                if kmer[pattern][-1]+k-kmer[pattern][0] <= L:
                    if pattern not in Patterns:
                        Patterns.append(pattern)
                del(kmer[pattern][0])
    return len(Patterns)
    
print(FindClumps(Genome,9,500,3))
#print(FindClumps('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA',5,50,4))
    



