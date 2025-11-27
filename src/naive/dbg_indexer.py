#!/usr/bin/env python3

import build
import query
import argparse
import pickle


def build_index(args):
    k = args.k
    file_list = args.i[0]
    out_name = args.o
    colored_graph = build.loop(file_list,k)
    pickle.dump(colored_graph,open(f"{out_name}.dumped","wb"))

def query_index(args):
    pass

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd')
    b = sub.add_parser('build')
    b.add_argument('-k', type=int, required=True, help= 'Size of each kmer')
    b.add_argument('-i', nargs='+', required=True,help='A list of fasta files, each representing a color')
    b.add_argument('-o', required=True, help= 'Name for the Cdbg you want to create')
    q = sub.add_parser('query')
    q.add_argument('--index', required=True)
    q.add_argument('--fasta', required=True, help='FASTA file with sequences to query')
    args = parser.parse_args()
    if args.cmd == 'build':
        build_index(args) # calling the main function for indexing
    elif args.cmd == 'query':
        query_index(args) # calling the main function for querying
    else:
        parser.print_help()
if __name__ == '__main__':
    main()