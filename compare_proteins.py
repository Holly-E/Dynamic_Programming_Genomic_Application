#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 27 09:04:26 2018

@author: hollyerickson

Compare the  Human Eyeless Protein and Fruitfly Eyeless Protein with PAX domain
"""
import matrix_functions as mf
import read_files as rf
import random
import numpy as np #use numpy to save dict
import matplotlib.pyplot as plt
import math

# Determine local alignment vals between Human and FruitFly Eyeless Proteins
local_align_matrix = mf.compute_alignment_matrix(rf.HumanEyelessProtein, rf.FruitflyEyelessProtein, rf.PAM50, False)
local_align_vals = mf.compute_local_alignment(rf.HumanEyelessProtein, rf.FruitflyEyelessProtein, rf.PAM50, local_align_matrix)

score =(local_align_vals[0]) #used in compute_z_score down below

def remove_dashes(str):
    str = str.replace("-","")
    return str

human_local = remove_dashes(local_align_vals[1])
fly_local = remove_dashes(local_align_vals[2])

# Determine global alignment between local alignments of Human & Fruitfly 
# Eyeless Proteins with the PAX Domain
human_gl_align_matrix = mf.compute_alignment_matrix(human_local, rf.ConsensusPAXDomain, rf.PAM50, True)
human_gl_align_vals = mf.compute_global_alignment(human_local, rf.ConsensusPAXDomain, rf.PAM50, human_gl_align_matrix) 
#print(human_gl_align_vals)

fly_gl_align_matrix = mf.compute_alignment_matrix(fly_local, rf.ConsensusPAXDomain, rf.PAM50, True)
fly_gl_align_vals = mf.compute_global_alignment(fly_local, rf.ConsensusPAXDomain, rf.PAM50, fly_gl_align_matrix) 
#print(fly_gl_align_vals)


def compute_percent_agree(seq_x, seq_y):
    total = len(seq_x)
    #print('len', total)
    
    if total != len(seq_y):
        return 'error in seq lengths'
    
    match = 0
    for index in range(total):
        if seq_x[index] == seq_y[index]:
            match += 1
    #print('match', match)       
    return match / total

#print (compute_percent_agree(human_gl_align_vals[1], human_gl_align_vals[2]))
#print (compute_percent_agree(fly_gl_align_vals[1], fly_gl_align_vals[2]))
  
    
def generate_null_distribution(seq_x,seq_y,scoring_matrix,num_trials):
    """
    Use statistical hypothesis testing to determine whether the local alignments
    computed are statistically significant.
    Returns a dictionary that represents an un-normalized distribution.
    """
    scoring_distribution = {}
   
    for trial in range(num_trials):
        rand_y = list(seq_y) # could also do [c for c in seq_y]
        random.shuffle(rand_y)
        rand_y = ''.join(rand_y) # different joining a str than javascript!
   
        # Determine local alignment score between seq_x and randomly ordered seq_y
        local_align_matrix = mf.compute_alignment_matrix(seq_x, rand_y, scoring_matrix, False)
        local_align_vals = mf.compute_local_alignment(seq_x, rand_y, scoring_matrix, local_align_matrix)
        score = local_align_vals[0]
        if score not in scoring_distribution.keys():
            scoring_distribution[score] = 1
        else:
            scoring_distribution[score] += 1
    
    return scoring_distribution
    
# scoring_dist = generate_null_distribution(rf.HumanEyelessProtein, rf.FruitflyEyelessProtein, rf.PAM50, 1000)
# np.save("d1.npy", scoring_dist)
scoring_dist_orig =np.load("d1.npy").item()
print(scoring_dist_orig)

def normalize_dist(dist, trials):
    """
    Normalize scoring distribution based off of #trials
    """
    y = [val / trials for key, val in dist.items()]
    
    return y

y = normalize_dist(scoring_dist_orig, 1000)
x = [key for key in scoring_dist_orig.keys()]

#plt.style.use('seaborn-whitegrid')
#fig = plt.figure(figsize = (9, 5))
#ax1 = fig.add_subplot(1, 1, 1)

#ax1.bar(x, y)

#ax1.set_title('Normalized Dist Over 1000 Trials', weight=600)
#ax1.set_xlabel('Scores')
#ax1.set_ylabel('Fraction of Total Trials')

#plt.savefig('Images/Normalized Dist', orientation='landscape')
#plt.show()

def mean(scoring_dist):
    scores = 0
    for key, val in scoring_dist.items():
        scores += key * val
    return scores / 1000

mean = mean(scoring_dist_orig)

def standard_deviation(scoring_dist, mean):
    scores_minus_mean = 0
    for key, val in scoring_dist.items():
        scores_minus_mean += (key - mean)**2 * val
    std_dev = math.sqrt(scores_minus_mean / 1000)
    return std_dev

std_dev = standard_deviation(scoring_dist_orig, mean)

def compute_z_score(score, mean, standard_dev):
    z_score = (score - mean) / standard_dev
    return z_score

z_score = compute_z_score(score, mean, std_dev) #score from line 21
#print(z_score) # z_score is number of standard deviations score is from mean



