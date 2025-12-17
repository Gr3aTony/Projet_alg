"""
Module query advanced

This module provides functions to query a colored De Bruijn graph
represented by unitigs using sequences from a FASTA file.
"""

from Bio import SeqIO


def query_similarity(cdbg: dict, seq: str, k: int, nbr_colors: int):
    """
    Compute similarity scores between a query sequence and a colored
    De Bruijn graph.

    For each k-mer extracted from the query sequence, the function checks
    whether it is contained in any unitig of the graph and updates
    similarity counts accordingly.

    Parameters:

    cdbg : dict
        Dictionary mapping unitigs to color identifiers.
    seq : str
        Query nucleotide sequence.
    k : int
        Length of the k-mers.
    nbr_colors : int
        Number of distinct colors in the graph.

    Returns:

    list of float
        Similarity scores for each color.
    """
    list_simili = [0] * nbr_colors

    for index in range(0, len(seq) - k):
        kmer = seq[index: index + k]

        for unitig in cdbg.keys():
            if kmer in unitig:
                for color in cdbg[unitig]:
                    list_simili[color] += 1
                break

    # Normalize similarity scores
    for i_list in range(nbr_colors):
        list_simili[i_list] = round(
            list_simili[i_list] / (len(seq) - k),
            4
        )

    return list_simili


def find_colors(cdbg):
    """
    Determine the number of distinct colors in a colored De Bruijn graph.

    Parameters:

    cdbg : dict
        Colored De Bruijn graph.

    Returns:

    int
        Number of distinct colors.
    """
    colors = set()

    for c in cdbg.values():
        for val in c:
            colors.add(val)

    return len(colors)


def query_compute(file_name: str, cdbg: dict, k: int):
    """
    Query a colored De Bruijn graph using sequences from a FASTA file.

    Parameters:

    file_name : str
        Path to a FASTA file containing query sequences.
    cdbg : dict
        Colored De Bruijn graph.
    k : int
        Length of the k-mers.

    Returns:

    str
        Formatted string containing sequence identifiers and similarity scores.
    """
    nbr_colors = find_colors(cdbg)
    out = ""

    for record in SeqIO.parse(file_name, "fasta"):
        seq = str(record.seq)
        similarity = query_similarity(cdbg, seq, k, nbr_colors)
        out += f"{record.id} {similarity}\n"

    return out
