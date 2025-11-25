import build
import query
import argparse
import pickle

#recup les args puis lancer build.loop(liste de fichier G ,taille kmer) stocker le dico puis creer l'objet
#parse et recup les info necessaire
def build_index(args):
    k = args.k
    file_list = args.i
    out_name = args.o
    colored_graph = build.loop(file_list,k)
    pickle.dump(colored_graph,open(f"{out_name}.dumped","wb"))

def query_index(args):
    pass

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd')
    b = sub.add_parser('build')
    b.add_argument('--k', type=int, required=True)
    b.add_argument('--fasta', nargs='+', required=True,
    help='FASTA files, one per color')
    b.add_argument('--out', required=True)
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