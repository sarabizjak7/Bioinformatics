import copy
import random 

import pickle
from os import path
from typing import Tuple, Generator, List


import sys
import threading


sys.setrecursionlimit(100000)
threading.stack_size(0x2000000)


##############################

def kmers(seq, k, stride = 1):
    """Generate a list of k-mers from a given string.

    Parameters
    ----------
    seq: str
        A string which will be decomposed into k-mers.
    k: int
    stride: int

    Returns
    -------
    List[str]

    Examples
    --------
    >>> kmers("mesenchyme", k=7, stride=1)
    ["mesench", "esenchy", "senchym", "enchyme"]

    >>> kmers("mesenchyme", k=7, stride=2)
    ["mesench", "senchym"]

    """
    list_of_kmers = []
    i = 0
    while i <= len(seq):
        if i + k <= len(seq):
            s = seq[i : i + k]
            list_of_kmers.append(s)
            i += stride
        else:
            break
    return list_of_kmers

##########################################################################

# Helper functions for Eulerian paths

def make_graph(seqs, k = None, stride = None):
     # Make a dict with nodes and connections in a graph
    graph = {}
    idx = 0
    for s in seqs:
        
        pref = s[0 : k - stride]
        suff = s[- (k - stride) :]

        if pref in graph:
            graph[pref].append((suff, s, idx))
        else:
            graph[pref] = [(suff, s, idx)]
        if suff not in graph:
            graph[suff] = []
        idx += 1

    return graph


def choose_starting_node(graph):
    counter = {edge: 0 for edge in graph}
    for edge1 in graph:
        for edge2 in graph:
            for conn in graph[edge2]:
                if edge1 in conn:
                    counter[edge1] += 1         
    for edge in counter:
        if counter[edge] == 0:
            starting_node = edge
            starting_seq = None
            break
        else:
            starting_node = list(graph.keys())[0]
            starting_seq = list(graph.values())[0][0][1]

    return starting_node, starting_seq


def possible_connections(v, graph):
    connections = graph[v]
    possible_conn = []
    for conn in connections:
        possible_conn.append((conn[0], conn[2]))
    return possible_conn


def connection(v, u, graph):
    for conn in graph[v]:
        if conn[0] == u:
            conn_val = (v, u, conn[2])
    return conn_val


def n_values(graph):
    count = 0
    for e in list(graph.values()):
        count += len(e)
    return count


def make_unique(l):
    unique_list = []
    for element in l:
        if element not in unique_list:
            unique_list.append(element)
        else:
            continue
    return unique_list


def make_tours(tours):
    T = []
    for tour in tours:
        t = []
        t.append(tour[0][0])
        t.append(tour[0][1])
        for conn in tour[1:]:
            t.append(conn[1])
        T.append(list(t))
    return make_unique(T)

 
def read_string(tours, graph, k, stride):
    all_string_tours = []
    for tour in tours:
        string_path = []
        for i in range(len(tour) - 1):
            node1 = tour[i]
            node2 = tour[i + 1]
            for conn in graph[node1]:
                if node2 in conn:
                    string_path.append(conn[1])
                    break     
        total_string = str(string_path[0])
        for s in string_path[1:]:
            total_string = total_string + s[k - stride:]
        all_string_tours.append(total_string)
    return all_string_tours

##############################

def assemble_genome (seqs, k = None, stride = None):
    """Perform genome assembly using the Eulerian path algorithm.

    The overall algorithm should follow the following general structure:
    1. For an input list of sequences, e.g. kmers, construct a DeBruijn graph as
       seen in the lectures.
    2. Find all possible Eulerian paths through the graph, i.e. all possible paths
       which visit each edge exactly once. Your paths should all start from a
       source node with in-degree zero. In case no such node exists, you may use
       the first sequence in the list of input sequences as your starting point.
    3. Decode your obtained paths into sequences, and return a list (or set) of
       unique genome assemblies as strings.

    Parameters
    ----------
    seqs: List[str]
        The list of strings. In this homework, this will always be a list of k-mers.

    k: int, optional
        This parameter is optional, and you may ignore it if you do not need it.
        But, this function will be called with the same `k` as `kmers`, which will
        be used to produce `seqs`. Therefore, make sure you do not remove this
        parameter.
    stride: int, optional
        This parameter is optional, and you may ignore it if you do not need it.
        But, this function will be called with the same `stride` as `kmers`, which will
        be used to produce `seqs`. Therefore, make sure you do not remove this
        parameter.

    Returns
    -------
    Set[str]
        A set (or list if you really want) of unique assemblies for the given `seqs`.

    """

    # Graph
    graph = make_graph(seqs, k, stride)

    # Copy graph
    G = copy.deepcopy(graph)

    # Choosing a node to start
    starting_node, starting_seq = choose_starting_node(graph)
   
    ##############################

    # Eulerian tour
    tours = []
    seen = set()
    
    def eulerian_tours(node, graph, tour):
        if n_values(G) > len(seen):
            neigh = possible_connections(node, graph)
            l = [v for v in neigh if (node, v[0], v[1]) not in seen]
      
            for next_node in l:
                e = (node, next_node[0], next_node[1])

                tour.append(e)
                seen.add(e)

                eulerian_tours(next_node[0], graph, tour)

                tour.pop()
                seen.remove(e)

        else:
            tours.append(list(tour))
    eulerian_tours(starting_node, G, [])

    # Write tours and read the string on the tour
    T = make_tours(tours)

    # upostevamo samo tiste, ki se zacnejo s prvim seq v seqs 
    if starting_seq != None:
        T_str = read_string(T, graph, k, stride)
        T_string = []
        for s in T_str:
            if s[:len(starting_seq)] == starting_seq:
                T_string.append(s)
    else:
        T_string = read_string(T, graph, k, stride)

    return T_string



