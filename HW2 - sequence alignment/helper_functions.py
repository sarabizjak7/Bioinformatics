import pickle
from os import path
from typing import Tuple, Generator, List

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy



####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

def load(organism_id: str) -> SeqRecord:
    """Load the NCBI record, use cached files if possible."""
    if not path.exists(path.join("data", f"{organism_id}.pkl.gz")):
        with Entrez.efetch(db="nucleotide", rettype="gb", id=organism_id) as handle:
            record = SeqIO.read(handle, "gb")
            with open(path.join("data", f"{organism_id}.pkl.gz"), "wb") as f:
                pickle.dump(record, f)
    else:
        with open(path.join("data", f"{organism_id}.pkl.gz"), "rb") as f:
            record = pickle.load(f)

    return record


####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################


def scoring(a, b):
    if a == "-" or b == "-":
        return -2
    if a != b:
        return -1
    else:
        return 2

    
def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    """
    # https://wilkelab.org/classes/SDS348/2019_spring/labs/lab13-solution.html?fbclid=IwAR10bQNBxLDvJ55odEXybijfzw6NtHK6hcaTHOM8apchxN8JBAmYbQwFuI


    indel = '-'
    seq1 = indel + seq1
    seq2 = indel + seq2

    n = len(seq1)
    m = len(seq2)
    
    # Score matrix

    score = numpy.zeros((n, m))

    # First row and first column

    score[0][0] = 0
    
    for i in range(n):
        score[i][0] = scoring_function(seq1[i], indel) * i
    for j in range(m):
        score[0][j] = scoring_function(seq2[j], indel) * j


    # Other values in the score matrix
    
    for i in range(1, n):
        for j in range(1, m):
            vertical = score[i - 1][j] + + scoring_function(seq1[i], indel)
            horizontal = score[i][j - 1] + scoring_function(indel, seq2[j])
            diagonal = score[i - 1][j - 1] + scoring_function(seq1[i], seq2[j])

            score[i][j] = max(vertical, horizontal, diagonal)


    align1 = ""
    align2 = ""

    i = n - 1
    j = m - 1

    while i > 0 and j > 0:
        current_score = score[i][j]
        vertical_score = score[i - 1][j]
        horizontal_score = score[i][j - 1]
        diagonal_score = score[i - 1][j - 1]

        if current_score == diagonal_score + scoring_function(seq1[i], seq2[j]):
            align1 += seq1[i]
            align2 += seq2[j]
            i -= 1
            j -= 1
        elif current_score == vertical_score + scoring_function(seq1[i], indel):
            align1 += seq1[i]
            align2 += indel
            i -= 1
        elif current_score == horizontal_score + scoring_function(indel, seq2[j]):
            align1 += indel
            align2 += seq2[j]
            j -= 1
            
    while j > 0:
        align1 += seq1[j-1]
        align2 += indel
        j -= 1
    while i > 0:
        align1 += indel
        align2 += seq2[i-1]
        i -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]
    optimal_score = score[n - 1][m - 1]

    return align1, align2, optimal_score



####################################################################################################################################################

def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    """

    indel = '-'

    seq1 = indel + seq1
    seq2 = indel + seq2

    n = len(seq1)
    m = len(seq2)
    
    # Score matrix

    score = numpy.zeros((n, m))

    # First row and first column

    score[0][0] = 0

    for i in range(1, n):
        score[i][0] = 0
    for j in range(1, m):
        score[0][j] = 0


    # Other values in the score matrix

    current_max = 0
    i_max = 0
    j_max = 0
    
    for i in range(1, n):
        for j in range(1, m):
            vertical = score[i - 1][j] + + scoring_function(seq1[i], indel)
            horizontal = score[i][j - 1] + scoring_function(indel, seq2[j])
            diagonal = score[i - 1][j - 1] + scoring_function(seq1[i], seq2[j])

            score[i][j] = max(vertical, horizontal, diagonal, 0)

            if score[i][j] >= current_max:
                current_max = score[i][j]
                i_max = i 
                j_max = j


    align1 = ""
    align2 = ""

    i = i_max

    j = j_max

    current_score = current_max

    while i > 0 and j > 0:
        #print("AAA")
        current_score = score[i][j]
        vertical_score = score[i - 1][j]
        horizontal_score = score[i][j - 1]
        diagonal_score = score[i - 1][j - 1]

        if current_score == diagonal_score + scoring_function(seq1[i], seq2[j]):
            #print("BBB")
            align1 += seq1[i]
            align2 += seq2[j]
            i -= 1
            j -= 1
        elif current_score == vertical_score + scoring_function(seq1[i], indel):
            #print("CCC")
            align1 += seq1[i]
            align2 += indel
            i -= 1
        elif current_score == horizontal_score + scoring_function(indel, seq2[j]):
            #print("DDD")
            align1 += indel
            align2 += seq2[j]
            j -= 1
            
        elif current_score == diagonal_score:
            #print("EEE")
            i -= 1
            j -= 1
        elif current_score == vertical_score:
            #print("FFF")
            i -= 1
        elif current_score == horizontal_score:
            #print("GGG")
            j -= 1
        elif current_score == 0:
            #print("break")
            break
        

    align1 = align1[::-1]
    align2 = align2[::-1]
    optimal_score = numpy.amax(score)

    return align1, align2, optimal_score
