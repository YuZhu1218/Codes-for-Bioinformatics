#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 00:36:48 2021

@author: yuzhu
"""


import networkx as nx
import pandas as pd
import numpy as np
import random
from collections import deque
import collections

# data=[]
# with open("dataset_203_7.txt","r") as f:
#     data = f.readlines()
# #print(data)
# for i in range(len(data)):
#     data[i] = data[i].strip('\n')
# Patterns = data[1::]
#print(Patterns)
    
    
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
    return G
    # adj_dict = {}   
    # for edge in G.edges:
    #     item1 = edge[0]
    #     item2 = edge[1]
    #     if item1 in adj_dict.keys():
    #         adj_dict[item1] += (',' + item2)
    #     else:
    #         adj_dict[item1] = item2
    # with open('6.txt', 'w') as f:
    #     for key in adj_dict.keys():
    #         f.write(key + ' -> ' + adj_dict[key] + '\n') 

def EulerianCycle(Graph):
    '''
    Code Challenge: Solve the Eulerian Cycle Problem.
    Input: The adjacency list of an Eulerian directed graph.
    Output: An Eulerian cycle in this graph.
    '''
    # G = nx.MultiDiGraph()
    # for i in range(len(Graph)):
    #     data = Graph[i].split(' -> ')
    #     G.add_node(data[0])
    #     out = data[1].split(',')
    #     for item in out:
    #         G.add_nodes_from(item)
    #         G.add_edge(data[0],item)
    G = Graph
    stack = deque()
    cycle = []
    current = random.choice(list(G.nodes()))
    stack.append(current)
    neighbors = list(G.neighbors(current))
    nex = random.choice(neighbors)
    G.remove_edge(current,nex)
    current = nex
    while len(stack) != 0:
        if list(G.neighbors(current)) == []:
            cycle.append(current)
            current = stack.pop()
        else:
            stack.append(current)
            neighbors = list(G.neighbors(current))
            nex = random.choice(neighbors)
            G.remove_edge(current,nex)
            current = nex
    cycle.append(cycle[0])
    cycle.reverse()
    return cycle

# Graph = ['0 -> 3',
# '1 -> 0',
# '2 -> 1,6',
# '3 -> 2',
# '4 -> 2',
# '5 -> 4',
# '6 -> 5,8',
# '7 -> 9',
# '8 -> 7',
# '9 -> 6']
#cycle = EulerianCycle(data)
#print(*cycle, sep='->')

def EulerianPath(Graph):
    '''
    Solve the Eulerian Path Problem.
    Input: The adjacency list of a directed graph that has an Eulerian path.
    Output: An Eulerian path in this graph.
    '''
    G = Graph
    # G = nx.MultiDiGraph()
    # for i in range(len(Graph)):
    #     data = Graph[i].split(' -> ')
    #     G.add_node(data[0])
    #     out = data[1].split(',')
    #     for item in out:
    #         G.add_nodes_from(item)
    #         G.add_edge(data[0],item)
    stack = deque()
    cycle = []
    start = []
    end = []
    for node in G.nodes():
        if G.in_degree(node) == (G.out_degree(node) + 1):
            end.append(node)
        elif (G.in_degree(node) + 1) == G.out_degree(node):
            start.append(node)
        elif G.in_degree(node) != G.out_degree(node):
            print('No Eulerian path')
            return 
    if (len(start) == 1) and (len(end) == 1):
        current = start[0]
    elif (start == []) and (end == []):
        return EulerianCycle(Graph)
    else: 
        print('No Eulerian path')
        return 
    stack.append(current)
    neighbors = list(G.neighbors(current))
    nex = random.choice(neighbors)
    G.remove_edge(current,nex)
    current = nex
    while len(stack) != 0:
        if list(G.neighbors(current)) == []:
            cycle.append(current)
            current = stack.pop()
        else:
            stack.append(current)
            neighbors = list(G.neighbors(current))
            nex = random.choice(neighbors)
            G.remove_edge(current,nex)
            current = nex
    cycle.append(start[0])
    cycle.reverse()
    return cycle

# Graph = ['0 -> 2',
# '1 -> 3',
# '2 -> 1',
# '3 -> 0,4',
# '6 -> 3,7',
# '7 -> 8',
# '8 -> 9',
# '9 -> 6']
# path = EulerianPath(data)
# print(*path, sep='->')
    
def StringReconstruction(k, Patterns):
    '''
    Solve the String Reconstruction Problem.
    Input: An integer k followed by a list of k-mers Patterns.
    Output: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may return any one.)

    '''
    dB = DeBrujinGraph(Patterns)
    path = EulerianCycle(dB)
    text = PathToGenome(path)
    return text

def IterateString(k):
    '''
    input: integer k
    output: list of k-mer binary integers
    '''
    Patterns = []
    for i in range(2**k):
        num = str(format(i, "b"))
        copy = ''
        for j in range(k-len(num)):
            copy += '0'
        copy += num
        Patterns.append(copy)
    return Patterns
        


# Patterns = ['CTTA',
# 'ACCA',
# 'TACC',
# 'GGCT',
# 'GCTT',
# 'TTAC',]
# Patterns = IterateString(8)
# print(StringReconstruction(8, Patterns))
    
def FindReadPairs(String,k,d):
    '''
    input: a stirng, k-mer with a gap d
    output: a dictionory, write on the txt
    '''
    length = len(String)
    dic = {}
    for i in range(length-2*k-d+1):
        dic[String[i:i+k]] = String[i+k+d:i+2*k+d]
    # od = collections.OrderedDict(sorted(dic.items()))
    return dic

od = FindReadPairs('TAATGCCATGGGATGTT',3,2)
for k,v in od.items():
    print('('+k+'|'+v+')')
# not work for repeated syntax #
    











    