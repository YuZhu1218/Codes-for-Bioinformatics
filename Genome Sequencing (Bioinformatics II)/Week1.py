#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 00:16:44 2021

@author: yuzhu
"""

import numpy as np
import pandas as pd
from pandas import DataFrame
import networkx as nx
import time


#Genome=remove(data[0])
#Genome=remove(data[1])
#print(Genome)

    
data=[]
with open("dataset_200_8.txt","r") as f:
    data = f.readlines()
#print(data)
for i in range(len(data)):
    data[i] = data[i].strip('\n')
#print(data)
#Genome=remove(data[0])
#Genome=remove(data[1])
#Genome = data[1]
#print(data)


def remove(string):
    return " ".join(string.split())


def FindGenome(data):
    Genome = []
    for i in range(1, len(data)):
        Genome.append(remove(data[i]))
    return Genome

#path = FindGenome(data)
#print(path)

def StringComposition(Genome,k):
    '''
    Generate the k-mer composition of a string.
    Input: An integer k and a string Text.
    Output: Compositionk(Text), where the k-mers are arranged in lexicographic order
    '''
    Strings = []
    for i in range(len(Genome)-k+1):
        Strings.append(Genome[i:i+k])
    return sorted(Strings)

# data = StringComposition(Genome,100)
# with open('1.txt', 'w') as f:
#     for item in data:
#         f.write(item)
#         f.write('\n')
        

def PathToGenome(path):
    '''
    String Spelled by a Genome Path Problem. Reconstruct a string from its genome path.
    Input: A sequence path of k-mers Pattern1, … ,Patternn such that the last k - 1 symbols of Patterni are equal to the first k-1 symbols
                of Patterni+1 for 1 ≤ i ≤ n-1.
    Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni (for 1 ≤ i ≤ n).
    '''
    String = path[0]
    for i in range(len(path)-1):
        String += path[i+1][-1]
    return String 

# path = ['ACCGA',
# 'CCGAA',
# 'CGAAG',
# 'GAAGC',
# 'AAGCT']
# print(PathToGenome(path))

# data = PathToGenome(data)
# with open('2.txt', 'w') as f:
#         f.write(data)
     

def OverlapGraph(Patterns):
    '''
    Input: A collection Patterns of k-mers.
    Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (You may return the nodes and and their edges  in any order)
    '''
    k = len(Patterns[0])
    dists = []
    for i in range(len(Patterns)):
        prefix = Patterns[i][0:k-1]
        suffix = Patterns[i][1:k]
        dists.append(DataFrame([1], index = [prefix], columns = [suffix]))
    adj_mat = pd.concat(dists).groupby(level=0).sum(min_count=1)
    adj_mat = adj_mat.fillna(0)
    return adj_mat

# adj_mat = OverlapGraph(data)
# with open('3.txt', 'w') as f:
#     for index in adj_mat.index:
#         adjacent = []
#         for column in adj_mat.columns:
#             if adj_mat.loc[index, column] == 1.0:
#                 adjacent.append(column)
#         f.write(str(index) + ' - > ' + str(','.join(adjacent)) + '\n')

def OverlapGraph2(Patterns):
    '''
    Input: A collection Patterns of k-mers.
    Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (You may return the nodes and and their edges  in any order)
    '''
    k = len(Patterns[0])
    G = nx.DiGraph()
    for i in range(len(Patterns)):
        G.add_node(Patterns[i])
    for node1 in G:
        for node2 in G:
            if node2[0:k-1] == node1[1:k]:
                G.add_edge(node1,node2)
    with open('4.txt', 'w') as f:
        for node in G:
            adjacency = G[node]
            if adjacency != {}:
                f.write(str(node) + ' -> ' + str(','.join(adjacency)) + '\n')

# start_time = time.time()         
# print(data)
# OverlapGraph2(data)
# print(f"time passed: {time.time() - start_time}")     

def DeBrujinGraphString(Text, k):
    '''
    De Bruijn Graph from a String Problem: Construct the de Bruijn graph of a string.
    Input: An integer k and a string Text.
    Output: DeBruijnk(Text).
    '''
    G = nx.MultiDiGraph()
    for i in range(len(Text)-k+1):
        node1 = Text[i:i+k-1]
        node2 = Text[i+1:i+k]
        G.add_node(node1)
        G.add_node(node2)
        G.add_edge(node1,node2)
    adj_dict = {}
    for edge in G.edges:
        item1 = edge[0]
        item2 = edge[1]
        if item1 in adj_dict.keys():
            adj_dict[item1] += (',' + item2)
        else:
            adj_dict[item1] = item2
    with open('5.txt', 'w') as f:
        for key in adj_dict.keys():
            f.write(key + ' -> ' + adj_dict[key] + '\n')

#Text = 'AAGATTCTCTAAGA'
# Text = 'TAATGCCATGGGATGTT'
# DeBrujinGraph(Text,2)
# DeBrujinGraph(Text,3)
# DeBrujinGraph(Text,4)
                
# G = nx.MultiDiGraph()
# G.add_node(1)
# G.add_node(2)
# G.add_edge(1,2)
# G.add_edge(1,2)
# print(G.edges())
            
def DeBrujinGraph(Patterns):
    '''
    DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.
    Input: A collection of k-mers Patterns.
    Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
    '''
    G = nx.MultiDiGraph()
    k = len(Patterns[0])
    for i in range(len(Patterns)):
        node1 = Patterns[i][0:k-1]
        node2 = Patterns[i][1:k]
        G.add_node(node1)
        G.add_node(node2)
        G.add_edge(node1,node2)
    adj_dict = {}   
    for edge in G.edges:
        item1 = edge[0]
        item2 = edge[1]
        if item1 in adj_dict.keys():
            adj_dict[item1] += (',' + item2)
        else:
            adj_dict[item1] = item2
    with open('6.txt', 'w') as f:
        for key in adj_dict.keys():
            f.write(key + ' -> ' + adj_dict[key] + '\n') 
# data = ['GAGG',
# 'CAGG',
# 'GGGG',
# 'GGGA',
# 'CAGG',
# 'AGGG',
# 'GGAG'] 
DeBrujinGraph(data)
    
    
    
    
    
    
    
    
    

        
