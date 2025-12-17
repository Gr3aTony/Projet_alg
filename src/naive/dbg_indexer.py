#!/usr/bin/env python3

import build
import query
import argparse
import pickle
import time


def build_index(args):
    k = args.k
    file_list = args.i[0]
    out_name = args.o

    debut_build = time.perf_counter()
    colored_graph = build.loop(file_list,k)
    fin_build = time.perf_counter()
    print(f"OUT TIME_BUILD: {fin_build - debut_build:.4f} second.")

    t0 = time.perf_counter()
    pickle.dump(colored_graph,open(f"{out_name}.dumped","wb"))
    t1 = time.perf_counter()
    print(f"OUT TIME_SERIALISATION: {t1 - t0:.4f} second.")

def query_index(args):
    q = args.q #query_file_name.fa

    t0 = time.perf_counter()
    cdbg = pickle.load(open( args.cdbg,"rb"))
    t1 = time.perf_counter()
    print(f"OUT TIME_DESERIALISATION {t1 - t0:.6f} second.")
    out_name = args.o# name of future output file

    debut_query = time.perf_counter()
    res = query.query_compute(q,cdbg)
    fin_query = time.perf_counter()
    print(f"OUT TIME_QUERY: {fin_query - debut_query:.4f} second.")

    with open(f"{out_name}.txt","w") as f:
        f.write(f"{res}")


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd')
    b = sub.add_parser('build')
    b.add_argument('-k', type=int, required=True, help= 'Size of each kmer')
    b.add_argument('-i', type=str, nargs='+', required=True,help='A list of fasta files, each representing a color')
    b.add_argument('-o', type=str, required=True, help= 'Name for the Cdbg you want to create')
    q = sub.add_parser('query')
    q.add_argument('-cdbg', required=True, help= 'A colored De Brujn Grpah of the Genome we want to compare the seq to')
    q.add_argument('-q', type=str, required=True, help='FASTA file with sequences to query') # might not be a dict type but a dumped obj
    q.add_argument('-o', type=str, required=True, help='Name of the desired output file')
    args = parser.parse_args()
    if args.cmd == 'build':
        build_index(args) # calling the main function for indexing
    elif args.cmd == 'query':
        query_index(args) # calling the main function for querying
    else:
        parser.print_help()
if __name__ == '__main__':
    #./dbg_indexer.py build -k 31 -i ../../G.txt -o trial1
    #./dbg_indexer.py query -cdbg trial1.dumped -q ../../First_set/query.fa -o ../../First_set/res_query_naive.txt
    main()