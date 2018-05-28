#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 26 16:47:24 2018

@author: hollyerickson

Function 1 builds a scoring matrix as a dictionary of dictionaries.
Function 2 computes an alignment matrix based off the scoring matrix.
Functions 3 & 4 use the alignment matrix returned by Function 2 to compute global and local alignments 
(respectively) of two sequences. 
"""

def build_scoring_matrix (alphabet,diag_score,off_diag_score,dash_score): 
    """
    Takes as input a set of characters and three scores and returns a dictionary of dictionaries
    whose entries are indexed by pairs of characters in alphabet plus -. 
    """
    scoring_matrix = {}
    
    scoring_matrix['-'] = {}
    scoring_matrix['-']['-'] = dash_score
    
    for index_i in alphabet:
        scoring_matrix[index_i] = {}
        scoring_matrix['-'][index_i] = dash_score
        scoring_matrix[index_i]['-'] = dash_score
        for index_j in alphabet:
            if index_i == index_j:
                scoring_matrix[index_i][index_j] = diag_score
            else:
                scoring_matrix[index_i][index_j] = off_diag_score
                
    return scoring_matrix

    
def compute_alignment_matrix(seq_x,seq_y,scoring_matrix,global_flag):
    """
    Takes seq_x andseq_y whose elements share a common alphabet with the scoring_matrix. 
    The function computes and returns the alignment matrix for seq_x and seq_y.  If global_flag is
    True, each entry of the alignment matrix is computed using the method described in Question 8 of the Homework. 
    If lobal_flag is False, each entry is computed using the method described in Question 12 of the Homework.
    """
    len_x = len(seq_x)
    len_y = len(seq_y)
    align_matrix = [[[] for _dummy_col in range(len_y + 1)] for _dummy_row in range(len_x + 1)] 
    align_matrix[0][0] = 0
    
    #populate matrix where seq_x and seq_y are only matched with dashes (row 0, col 0)
    for index_i in range(1, len_x + 1):
        score = align_matrix[index_i - 1][0] + scoring_matrix[seq_x[index_i - 1]]['-']
        # assign score of 0 if computing local alignment with negative score
        if score >= 0 or global_flag == True:
            align_matrix[index_i][0] = score
        else:
            align_matrix[index_i][0] = 0
            
    for index_j in range(1, len_y + 1):
        score = align_matrix[0][index_j - 1] + scoring_matrix['-'][seq_y[index_j - 1]]
        if score >= 0 or global_flag == True:
            align_matrix[0][index_j] = score
        else:
            align_matrix[0][index_j] = 0
    
    for index_i in range(1, len_x + 1):
        for index_j in range(1, len_y + 1):
            diag = align_matrix[index_i - 1][index_j - 1] + scoring_matrix[seq_x[index_i - 1]][seq_y[index_j - 1]]
            top = align_matrix[index_i - 1][index_j] + scoring_matrix[seq_x[index_i - 1]]['-']
            left = align_matrix[index_i][index_j - 1] + scoring_matrix['-'][seq_y[index_j - 1]]
            score = max(diag, top, left)
            
            if score >= 0 or global_flag == True:
                align_matrix[index_i][index_j] = score
            else:
                align_matrix[index_i][index_j] = 0
    return align_matrix


def compute_global_alignment(seq_x,seq_y,scoring_matrix,align_matrix): 
    """
    Takes seq_x and seq_y whose elements share a common alphabet with the scoring matrix.
    Computes a global alignment of the sequences using the global alignment matrix. 
    Returns a tuple of the form (score,align_x,align_y). Note that align_x and align_y should have the 
    same length and may include the padding character -.
    """
    index_i = len(seq_x)
    index_j = len(seq_y)
    align_x = ''
    align_y = ''
    score = align_matrix[index_i][index_j]
    
    while index_i != 0 and index_j != 0:
        if align_matrix[index_i][index_j] == align_matrix[index_i - 1][index_j - 1] + scoring_matrix[seq_x[index_i - 1]][seq_y[index_j - 1]]:
             # if current position came from diagonal 
            align_x = seq_x[index_i -1] + align_x
            align_y = seq_y[index_j - 1] + align_y
            index_i -= 1
            index_j -= 1
        elif align_matrix[index_i][index_j] == align_matrix[index_i - 1][index_j] + scoring_matrix[seq_x[index_i - 1]]['-']:
            # if current pos came from top
            align_x = seq_x[index_i -1] + align_x
            align_y = '-' + align_y
            index_i -= 1
        else:
            # if current pos came from left
            align_x = '-' + align_x
            align_y = seq_y[index_j - 1] + align_y
            index_j -= 1
    
    # consider row 0 and col 0 
    while index_i != 0:
        align_x = seq_x[index_i -1] + align_x
        align_y = '-' + align_y
        index_i -= 1
        
    while index_j != 0:
        align_x = '-' + align_x
        align_y = seq_y[index_j - 1] + align_y
        index_j -= 1
        
    return (score, align_x, align_y)

def compute_local_alignment(seq_x,seq_y,scoring_matrix,align_matrix):
    """
    Similar to compute_global_alignment except computes an optimal local alignment 
    starting at the maximum entry of the local alignment matrix and working backwards.
    The align sequences will stop when the alignment matrix reaches a score of zero.
    """
    #get starting position (max entry of local alignment)
    score = -1
    start_pos = [0, 0]
    for pos_x in range(0, len(seq_x) + 1):
        for pos_y in range(0, len(seq_y) + 1):
            if align_matrix[pos_x][pos_y] > score:
                score = align_matrix[pos_x][pos_y]
                start_pos = [pos_x, pos_y]
    
    index_i = start_pos[0]
    index_j = start_pos[1]
    align_x = ''
    align_y = ''
             
    while index_i != 0 and index_j != 0 and align_matrix[index_i][index_j] != 0:
        if align_matrix[index_i][index_j] == align_matrix[index_i - 1][index_j - 1] + scoring_matrix[seq_x[index_i - 1]][seq_y[index_j - 1]]:
            # if current position came from diagonal 
            align_x = seq_x[index_i -1] + align_x
            align_y = seq_y[index_j - 1] + align_y
            index_i -= 1
            index_j -= 1
        elif align_matrix[index_i][index_j] == align_matrix[index_i - 1][index_j] + scoring_matrix[seq_x[index_i - 1]]['-']:
            # if current pos came from top
            align_x = seq_x[index_i -1] + align_x
            align_y = '-' + align_y
            index_i -= 1
        else:
            # if current pos came from left
            align_x = '-' + align_x
            align_y = seq_y[index_j - 1] + align_y
            index_j -= 1
    
    # consider row 0 and col 0 
    while index_i != 0 and align_matrix[index_i][index_j] != 0:
        align_x = seq_x[index_i -1] + align_x
        align_y = '-' + align_y
        index_i -= 1
        
    while index_j != 0 and align_matrix[index_i][index_j] != 0:
        align_x = '-' + align_x
        align_y = seq_y[index_j - 1] + align_y
        index_j -= 1
        
    return (score, align_x, align_y)

    
#score_matrix = build_scoring_matrix('abcd', 10, 4, -6)
#align_matrix = compute_alignment_matrix('ab','abd', score_matrix, False)
#print(align_matrix)
#print(compute_local_alignment('ab','abd', score_matrix, align_matrix))
#print (compute_local_alignment('ACTACT', 'AGCTA', {'A': {'A': 2, 'C': 1, '-': 0, 'T': 1, 'G': 1}, 'C': {'A': 1, 'C': 2, '-': 0, 'T': 1, 'G': 1}, '-': {'A': 0, 'C': 0, '-': 0, 'T': 0, 'G': 0}, 'T': {'A': 1, 'C': 1, '-': 0, 'T': 2, 'G': 1}, 'G': {'A': 1, 'C': 1, '-': 0, 'T': 1, 'G': 2}}, [[0, 0, 0, 0, 0, 0], [0, 2, 2, 2, 2, 2], [0, 2, 3, 4, 4, 4], [0, 2, 3, 4, 6, 6], [0, 2, 3, 4, 6, 8], [0, 2, 3, 5, 6, 8], [0, 2, 3, 5, 7, 8]]))