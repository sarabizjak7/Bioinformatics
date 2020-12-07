from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from os import path
import pickle

import numpy as np

# Source : https://github.com/lex8erna/UPGMApy/blob/master/UPGMA.py

def UPGMA(distances):
    """
    Unweighted pair group method with arithmetic mean (UPGMA) agglomerative clustering.

    Parameters
    ----------
    distances: np.ndarray
        A two dimensional, square, symmetric matrix containing distances between data
        points. The diagonal is zeros.

    Returns
    -------
    np.ndarray
        The linkage matrix, as specified in scipy. Briefly, this should be a 2d matrix
        each row containing 4 elements. The first and second element should denote the
        cluster IDs being merged, the third element should be the distance, and the
        fourth element should be the number of elements within this new cluster. Any
        new cluster should be assigned an incrementing ID, e.g. after the first step
        where the first two points are merged into a cluster, the new cluster ID should
        be N, then N+1, N+2, ... in subsequent steps.
    """

    n = len(distances)

    # set the output matrix
    M_output = np.zeros((n - 1, 4))

    cluster_IDs = [[i] for i in range(n)]

    # (n - 1) times : we search for min dist in distance matrix and then join them into one new ID
    for i in range(n - 1):

        ##################
        # min distance in input matrix and its index
        min_distance = float("inf")
        ind_x = -1
        ind_y = -1

        for x in range(len(distances)):
            for y in range(x + 1, len(distances[x])):
                if distances[x][y] < min_distance:
                    min_distance = distances[x][y]
                    ind_x = x
                    ind_y = y
        ##################

        # combine clusters x and y where distance is minimal
        # we save the new cluster into cluster_IDs as a list of combined clusters, so we can track the length
        new_cluster = cluster_IDs[ind_x] + cluster_IDs[ind_y]
        cluster_IDs.append(new_cluster)

        # update the output matrix
        M_output[i][0] = ind_x
        M_output[i][1] = ind_y
        M_output[i][2] = min_distance
        M_output[i][3] = len(new_cluster)

        # expand the distance matrix
        distances = np.append(distances, np.atleast_2d(np.zeros(n + i)), axis=0)
        #print("D1: " + str(distances))
        distances = np.hstack((distances, np.atleast_2d(np.zeros(n + i + 1)).T))
        #print("D2: " + str(distances))

        # we add the new cluster to the matrix by expanding it
        # we are on i-th iteration, so the matrix is of len (n + i)

        for j in range(n + i + 1):
            # distances with all other clusters
            for k in range(j + 1, n + i + 1):
                if k == n + i and j != ind_x and j != ind_y:
                    # calculations of new distances:
                    d = (distances[j][ind_x] * len(cluster_IDs[ind_x]) + distances[j][ind_y] * len(cluster_IDs[ind_y])) / (len(cluster_IDs[ind_x]) + len(cluster_IDs[ind_y]))
                    # update distance matrix
                    distances[j][k] = d
                    distances[k][j] = d

        # "delete" rows and columns that we combined
        distances[:, ind_x] = np.nan
        distances[:, ind_y] = np.nan
        distances[ind_x, :] = np.nan
        distances[ind_y, :] = np.nan
        
    return M_output


def jukes_cantor(p: float) -> float:
    """The Jukes-Cantor correction for estimating genetic distances.

    Parameters
    ----------
    p: float
        The proportional distance, i.e. the number of of mismatching symbols (Hamming
        distance) divided by the total sequence length.

    Returns
    -------
    float
        The corrected genetic distance.

    """
    return -(3 / 4) * np.log(1 - 4 / 3 * p)




##########################################

def find_orfs(sequence, start_codons, stop_codons):
    """Find possible ORF candidates in a single reading frame.
    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]
    Returns
    -------
    List[Tuple[int, int]]
        tuples of form (start_loc, stop_loc)
    """
    
    ORF_candidates = []
    start_pos = []
    stop_pos = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i: i + 3]
        if codon in start_codons:
            start_pos.append(i)
        if codon in stop_codons:
            stop_pos.append(i + 3) # Stop codon is included
            
    #print("START:" + str(start_pos))
    #print("STOP:" + str(stop_pos))

    # Append cantidates for ORFS in a prepared list
    stop_p = 0
    for start in start_pos:
        if start >= stop_p :
            for stop in stop_pos:   
                if stop > start:
                    # Update the stop position
                    stop_p = stop
                    ORF_candidates.append((start, stop))
                    break
            
    return ORF_candidates
        

####################################################################################################################################################

    
def find_all_orfs(sequence, start_codons, stop_codons):
    """Find ALL the possible ORF candidates in the sequence using all six
    reading frames.
    Parameters
    ----------
    sequence: Seq
    start_codons: List[str]
    stop_codons: List[str]
    Returns
    -------
    List[Tuple[int, int, int]]
        tuples of form (strand, start_loc, stop_loc). Strand should be either 1
        for reference strand and -1 for reverse complement.
    """
    all_ORF_candidates = []

    # Reference strand (starting at pos 0 (sequence), starting at pos 1, starting at pos 2)
    # Find ORFs for each one of three sequences
    strand0 = find_orfs(sequence, start_codons, stop_codons)
    strand1 = find_orfs(sequence[1:], start_codons, stop_codons)
    strand2 = find_orfs(sequence[2:], start_codons, stop_codons)
    #print(strand0)
    #print(strand1)
    #print(strand2)
    
    # Add ORFs to a list with all candidates -- tuples with 1 
    [all_ORF_candidates.append((1, start0, stop0)) for (start0, stop0) in strand0]
    [all_ORF_candidates.append((1, start1 + 1, stop1 + 1)) for (start1, stop1) in strand1]
    [all_ORF_candidates.append((1, start2 + 2, stop2 + 2)) for (start2, stop2) in strand2]

    #########################

    rev_seq = sequence.reverse_complement()

    # Reversed stand (starting at pos 0 (sequence), starting at pos 1, starting at pos 2)
    # Find ORFs for each one of three reversed sequences ([::-1])
    rev_strand0 = find_orfs(rev_seq, start_codons, stop_codons)
    rev_strand1 = find_orfs(rev_seq[1:], start_codons, stop_codons)
    rev_strand2 = find_orfs(rev_seq[2:], start_codons, stop_codons)
    #print(rev_strain0)
    #print(rev_strain1)
    #print(rev_strain2)
    
    # Add ORFs to a list with all candidates -- tuples with -1
    # Indexes as in reference strand -> len(seq) - index_of_reversed_strand 
    [all_ORF_candidates.append((-1, len(sequence) - stop0, len(sequence) - start0)) for (start0, stop0) in rev_strand0]
    [all_ORF_candidates.append((-1, len(sequence) - (stop1 + 1), len(sequence) - (start1 + 1))) for (start1, stop1) in rev_strand1]
    [all_ORF_candidates.append((-1, len(sequence) - (stop2 + 2), len(sequence) - (start2 + 2))) for (start2, stop2) in rev_strand2]

    return all_ORF_candidates