"""
Module build naive

This module provides functions to build a colored De Bruijn graph from
multiple FASTA files.

Each extracted k-mer is associated with one or more colors, representing
the origin file in which the k-mer appears.
"""

from Bio import SeqIO


def dbg_builder(fasta_file: str, color: int, k: int, dbg: dict):
    """
    Update a colored De Bruijn graph using a single FASTA file.

    This function parses a FASTA file, extracts all possible k-mers from
    its sequences, and inserts them into a colored De Bruijn graph.
    Each k-mer is associated with a color corresponding to the input file.

    Parameters :
  
    fasta_file : str
        Path to a FASTA file containing one or more sequences.
    color : int
        Integer identifier representing the color associated with the FASTA file.
    k : int
        Length of the k-mers to extract.
    dbg : dict
        Dictionary representing the colored De Bruijn graph, where keys
        are k-mers (str) and values are lists of color identifiers.

        
    Returns :

    dict
        Updated colored De Bruijn graph including k-mers from the FASTA file.
    """
    # Iterate over all sequences contained in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Explicit conversion of the sequence to a string
        seq = str(record.seq)

        # Sliding window extraction of k-mers along the sequence
        for i in range(0, len(seq) - k):
            kmer = seq[i: i + k]

            # Add the k-mer if it does not yet exist in the graph
            if kmer not in dbg.keys():
                dbg[kmer] = [color]

            # Add the color if the k-mer exists but is not yet associated with it
            elif color not in dbg[kmer]:
                dbg[kmer].append(color)

    return dbg


def loop(file_list: str, k: int):
    """
    Build a colored De Bruijn graph from a list of FASTA files.

    This function reads a text file containing paths to multiple FASTA files
    (one per line) and iteratively updates a shared De Bruijn graph by calling
    `dbg_builder` on each file.

    Parameters :

    file_list : str
        Path to a text file containing one FASTA file path per line.
    k : int
        Length of the k-mers to extract.

        
    Returns :
   
    dict
        Colored De Bruijn graph built from all input FASTA files.
    """
    dbg = {}
    i = 0

    # Open the file listing the FASTA file paths
    with open(file_list, "r") as fl:
        for file in fl:
            # Remove trailing whitespace and newline characters
            file = file.strip()

            # Update the graph with k-mers from the current file
            dbg = dbg_builder(file, i, k, dbg)

            # Increment the color identifier
            i += 1

    return dbg
