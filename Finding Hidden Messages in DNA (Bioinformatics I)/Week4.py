#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 21:24:31 2021

@author: yuzhu
"""

import random
import Week3 as wk

data=[]
with open("datase7.txt","r") as f:
    data = f.readlines()
for i in range(len(data)):
    data[i] = data[i].strip('\n')

def RandomizedMotifSearch(Dna, k, t):
    Motifs = []
    for i in range(len(Dna)):
        n = len(Dna[i])-k
        num = random.randint(0,n)
        Motifs.append(Dna[i][num:num+k])
    BestMotifs = Motifs[:]
    while True:
        Profile = wk.formProfile(BestMotifs)
        Motifs = []
        for j in range(len(Dna)):
            pattern = wk.ProfileMostProbableKmer(Dna[j], k, Profile)
            Motifs.append(pattern)
        if wk.CountEntropy(Motifs) < wk.CountEntropy(BestMotifs):
            BestMotifs = Motifs[:]
        else:
            return BestMotifs
    
def ThousandRandomizedMotifSearch(Dna, k, t):
    for i in range(1000):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        if i == 0:
            BestMotifs = Motifs[:]
        elif wk.CountEntropy(Motifs) < wk.CountEntropy(BestMotifs):
           BestMotifs = Motifs[:]
    return BestMotifs

# Dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
# 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
# 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
# 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
# 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']

# answer = ['TCTCGGGG',
# 'CCAAGGTG',
# 'TACAGGCG',
# 'TTCAGGTG',
# 'TCCACGTG']

Dna = data[1::]
print(Dna)
# result = ThousandRandomizedMotifSearch(Dna, 15, 20)
# print(*result, sep='\n')
# #print(wk.CountEntropy(result))
# #print(wk.CountEntropy(answer))


def GibbsSampler(Dna, k, t, N):
    Motifs = []
    for i in range(len(Dna)):
        n = len(Dna[i])-k
        num = random.randint(0,n)
        Motifs.append(Dna[i][num:num+k])
    BestMotifs = Motifs[:]
    for j in range(N):
        MotifsCopy = Motifs[:]
        i = random.randint(0,t-1)
        del(MotifsCopy[i])
        Profile = wk.formProfile(MotifsCopy)
        Prob = []
        for n in range(len(Dna[i])-k+1):
            Pattern = Dna[i][n:n+k]
            Number = 1
            for m in range(int(k)):
                base = Pattern[m]
                Number *= Profile[base][m]
            Prob.append(Number)
        MotifiNum = random.choices(range(len(Dna[i])-k+1),weights = Prob, k=1)
        Motifi = Dna[i][MotifiNum[0]:MotifiNum[0]+k]
        Motifs[i] = Motifi
        if wk.CountEntropy(Motifs) < wk.CountEntropy(BestMotifs):
            BestMotifs = Motifs[:]
    return BestMotifs

# Dna = ['CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA',
# 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
# 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
# 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
# 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']

result = GibbsSampler(Dna, 15, 20, 2000)
print(*result, sep=' ')
print(wk.CountEntropy(result))