#!/usr/bin/env python3
"""
Module dbg_indexer advanced

This module implements a command-line interface for building and querying
a colored De Bruijn graph.

It provides two main subcommands:
- buil`: constructs a colored De Bruijn graph from input FASTA files
  and serializes it to disk;
- query: loads a serialized graph and queries it using FASTA sequences.
"""

import build
import query
import argparse
import pickle
import time


def build_index(args):
    """
    Build and serialize a colored De Bruijn graph.

    This function constructs a colored De Bruijn graph using the input
    FASTA files specified in the command-line arguments, measures the
    construction time, and serializes the resulting graph to disk.

    Parameters:

    args : argparse.Namespace
        Parsed command-line arguments containing:
        - k: size of the k-mers,
        - i: list containing the path to a file listing FASTA files,
        - o: output file prefix for the serialized graph.
    """
    k = args.k
    file_list = args.i[0]
    out_name = args.o

    # Measure the time required to build the colored De Bruijn graph
    debut_build = time.perf_counter()
    colored_graph = build.loop(file_list, k)
    fin_build = time.perf_counter()
    print(f"OUT TIME_BUILD: {fin_build - debut_build:.6f} second.")

    # Serialize the graph to disk and measure serialization time
    t0 = time.perf_counter()
    pickle.dump(colored_graph, open(f"{out_name}.dumped", "wb"))
    t1 = time.perf_counter()
    print(f"OUT TIME_SERIALISATION: {t1 - t0:.6f} second.")


def query_index(args):
    """
    Query a serialized colored De Bruijn graph.

    This function loads a previously serialized colored De Bruijn graph,
    measures deserialization time, performs a query using FASTA sequences,
    and writes the results to an output file.

    Parameters:
    
    args : argparse.Namespace
        Parsed command-line arguments containing:
        - cdbg: path to the serialized colored De Bruijn graph,
        - q: FASTA file containing query sequences,
        - o: output file name.
    """
    q = args.q  # FASTA file containing query sequences

    # Load the serialized graph and associated k-mer bag
    t0 = time.perf_counter()
    cdbg, kmer_bag = pickle.load(open(args.cdbg, "rb"))
    t1 = time.perf_counter()
    print(f"OUT TIME_DESERIALISATION: {t1 - t0:.6f} second.")

    # Infer k-mer length from the stored k-mer bag
    k = len(next(iter(kmer_bag)))
    out_name = args.o

    # Perform the query and measure execution time
    debut_query = time.perf_counter()
    res = query.query_compute(q, cdbg, k)
    fin_query = time.perf_counter()
    print(f"OUT TIME_QUERY: {fin_query - debut_query:.6f} second.")

    # Write query results to the output file
    with open(f"{out_name}.txt", "w") as f:
        f.write(f"{res}")


def main():
    """
    Parse command-line arguments and dispatch subcommands.

    This function defines the command-line interface, parses user input,
    and calls the appropriate function depending on the selected subcommand.
    """
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd')

    # Subcommand for building a colored De Bruijn graph
    b = sub.add_parser('build')
    b.add_argument('-k', type=int, required=True, help='Size of each k-mer')
    b.add_argument(
        '-i',
        type=str,
        nargs='+',
        required=True,
        help='A list of FASTA files, each representing a color'
    )
    b.add_argument(
        '-o',
        type=str,
        required=True,
        help='Name for the colored De Bruijn graph to create'
    )

    # Subcommand for querying a colored De Bruijn graph
    q = sub.add_parser('query')
    q.add_argument(
        '-cdbg',
        required=True,
        help='Serialized colored De Bruijn graph to query'
    )
    q.add_argument(
        '-q',
        type=str,
        required=True,
        help='FASTA file containing query sequences'
    )
    q.add_argument(
        '-o',
        type=str,
        required=True,
        help='Name of the desired output file'
    )

    args = parser.parse_args()

    # Dispatch according to the selected subcommand
    if args.cmd == 'build':
        build_index(args)
    elif args.cmd == 'query':
        query_index(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
