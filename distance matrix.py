#!/usr/bin/env python
# coding: utf-8

# In[209]:


from numpy import where
import numpy as np
import pandas as pd
from numpy import fill_diagonal
from pandas import DataFrame
from itertools import combinations
import glob
import sys
from Bio import SeqIO
import json


'''
this code will explain how to compute the distance matrix for multiple sequence aligen data.
'''

def pairwise_distance(first_seq, second_seq):
    '''
    compute the pairwise distances between two sequences
    '''
    total_distance = 0
    
    if len(first_seq) > len(second_seq):
        total_distance =  len(first_seq) - len(second_seq)
    elif len(second_seq) > len(first_seq):
        total_distance = len(second_seq) - len(first_seq)
    
    
    for char_1, char_2 in zip(first_seq, second_seq):
        if char_1 != char_2:
            total_distance += 1
    return total_distance


def distance_matrix(dist_dict):
    '''
   Create a blank dataframe and label its column and index headings
   with all of the creatures.
    '''
    dist_matrix = DataFrame(index=dist_dict.keys(), columns=dist_dict.keys())
    fill_diagonal(dist_matrix.values, 0)
    return dist_matrix


def generate_distance_matrix(dist_matrix, pairwise_dict, start, counter, end, sequence_dict, sequence_df):
    '''
    Utilize vectorized assignments to fill in the complete distance matrix.
    Imagine allocating lists for each organism that are at right angles and have a vertex of zero.
    Setting the start position to the end position and adding a decrementing counter value to the
    end value update the start and end positions for the right angles.
    '''
    for i, j in zip(range(len(sequence_dict.keys()) - 1), range(1, len(sequence_dict.keys()))):
        sequence_df.iloc[i, j:] = list(pairwise_dict.values())[start:end]
        sequence_df.iloc[j:, i] = sequence_df.iloc[i, j:]
        start = end
        end += counter
        counter -= 1
    return dist_matrix



# In[210]:


def generate_distance_matrix_for_genotype(data):
    '''
     input : this block of code will take the JSON dataset which contains the sample name and their genotype data 
     which is present in array format 

     output : this block would generate the distance matrix of multiple sequence 

    '''
    df = pd.DataFrame(data["numbers_of_alternate_alleles"], index=data["samples_selected_mapped"])
    df = df.astype(str)
    columns=df.index
    sequence_dict = {}

    for x in range(len(columns)):

        sequence_dict[columns[x]] = str(''.join(df.iloc[x]))

    sequence_df = distance_matrix(sequence_dict)
    Seq_pairwise_diff = {
        key: pairwise_distance(sequence_dict[key[0]], sequence_dict[key[1]])
        for key in combinations(sequence_dict, 2)
    }

    Seq_df = generate_distance_matrix(sequence_df, Seq_pairwise_diff, 0,
                            len(sequence_dict.keys()) - 2,
                            len(sequence_dict.keys()) - 1,sequence_dict,sequence_df)
    
    return Seq_df


# In[211]:


f = open('data/sample.json')
data = json.load(f)
Seq_df = generate_distance_matrix_for_genotype(data)


# In[212]:


Seq_df


# In[213]:


def generate_distance_matrix_for_sequence(all_fasta_file):
    '''
        input : this block of code will take the fattest dataset which contains the sample name and their sequence data 
            which is present in array format

        output : this block would generate the distance matrix of multiple sequence 

    '''

    sequence_list = [] # To keep order of sequence
    sequence_dict = {}
    for record in SeqIO.parse(open(all_fasta_file, "r"), "fasta"):
        sequence_list.append(record.id)
        sequence_dict[record.id] = str(record.seq)

    sequence_df = distance_matrix(sequence_dict)
    Seq_pairwise_diff = {
        key: pairwise_dist(sequence_dict[key[0]], sequence_dict[key[1]])
        for key in combinations(sequence_dict, 2)
    }

    Seq_df = generate_distance_matrix(sequence_df, Seq_pairwise_diff, 0,
                                len(sequence_dict.keys()) - 2,
                                len(sequence_dict.keys()) - 1, sequence_dict,sequence_df)
    
    return Seq_df


# In[214]:


all_fasta_file = "data/vcf_to_seq.fasta"
Seq_df = generate_distance_matrix_for_sequence(all_fasta_file)


# In[215]:


Seq_df

