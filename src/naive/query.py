"""
Module query naive

This module provides functions to query a colored De Bruijn graph using
sequences from a FASTA file.

It computes similarity scores between query sequences and the genomes
represented in the colored De Bruijn graph.
"""

from Bio import SeqIO


def query_similarity(cdbg: dict, seq: str, k: int, nbr_colors: int):
    """
    Compute similarity scores between a sequence and a colored De Bruijn graph.

    This function extracts all k-mers from the input sequence and counts
    how many of them are present in the colored De Bruijn graph for each
    color. Similarity scores are expressed as proportions.


    Parameters:

    cdbg : dict
        Colored De Bruijn graph represented as a dictionary where keys are
        k-mers and values are lists of color identifiers.
    seq : str
        Nucleotide sequence to query against the graph.
    k : int
        Length of the k-mers used in the graph.
    nbr_colors : int
        Total number of distinct colors in the graph.


    Returns:
    list of float
        List containing similarity scores for each color.
    """
    # Initialize similarity counters for each color
    list_simili = [0] * nbr_colors

    # Slide a window of size k along the sequence
    for index in range(0, len(seq) - k):
        kmer = seq[index: index + k]

        # Check whether the k-mer exists in the colored De Bruijn graph
        if kmer in cdbg:
            # Increment the similarity counter for each associated color
            for color in cdbg[kmer]:
                list_simili[color] += 1

    # Normalize counts to obtain similarity proportions
    for i_list in range(nbr_colors):
        list_simili[i_list] = round(
            list_simili[i_list] / (len(seq) - k),
            4
        )

    return list_simili


def find_colors(cdbg):
    """
    Determine the number of distinct colors in a colored De Bruijn graph.

    This function scans the values of the colored De Bruijn graph and
    collects all unique color identifiers.


    Parameters:

    cdbg : dict
        Colored De Bruijn graph represented as a dictionary.


    Returns:

    int
        Number of distinct colors found in the graph.
    """
    colors = set()

    # Collect all unique color identifiers from the graph
    for c in cdbg.values():
        for val in c:
            colors.add(val)

    return len(colors)


def query_compute(file_name: str, cdbg: dict):
    """
    Query a colored De Bruijn graph using sequences from a FASTA file.

    This function parses a FASTA file, computes similarity scores between
    each sequence and the colored De Bruijn graph, and formats the results
    as a text output.
    

    Parameters:

    file_name : str
        Path to a FASTA file containing query sequences.
    cdbg : dict
        Colored De Bruijn graph to query.

    Returns :
    str
        Formatted string containing sequence identifiers and similarity
        scores for each color.
    """
    # Infer the k-mer size from the keys of the graph
    k = len(next(iter(cdbg)))

    # Determine the number of distinct colors in the graph
    nbr_colors = find_colors(cdbg)

    out = ""

    # Iterate over query sequences from the FASTA file
    for record in SeqIO.parse(file_name, "fasta"):
        seq = str(record.seq)

        # Compute similarity scores for the current sequence
        similarity = query_similarity(cdbg, seq, k, nbr_colors)

        # Append formatted results to the output string
        out += f"{record.id} {similarity}\n"

    return out
