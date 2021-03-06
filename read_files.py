#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 27 01:21:57 2018

@author: hollyerickson

For the genomics part, load several protein sequences and an appropriate scoring matrix. 
For the spelling correction part of the Application, load a provided word list.
"""
#urllib2 not available in Python3
from urllib.request import urlopen


# URLs for data files
PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"



###############################################
# provided code

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = urlopen(filename)
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
        vals = line.split()
        xkey = vals.pop(0).decode("utf-8") #add if opening in py3 to remove bytes prefix 
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            ykey = ykey.decode("utf-8") # must be decoded at each element within a list, not the full list!!
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict


def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = urlopen(filename)
    protein_seq = protein_file.read().decode("utf-8") #add if opening in py3
    protein_seq = protein_seq.rstrip()
    return protein_seq


def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = urlopen(filename)
    
    # read in files as string
    words = word_file.read().decode("utf-8") #add if opening in py3
    
    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print ("Loaded a dictionary with", len(word_list), "words")
    return word_list

PAM50 = read_scoring_matrix(PAM50_URL) 

HumanEyelessProtein = read_protein(HUMAN_EYELESS_URL)
FruitflyEyelessProtein  = read_protein(FRUITFLY_EYELESS_URL)
ConsensusPAXDomain = read_protein(CONSENSUS_PAX_URL)

word_list = read_words(WORD_LIST_URL)