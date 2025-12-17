from Bio import SeqIO
def query_similarity(cdbg : dict, seq : str , k : int, nbr_colors : int):
    """
    Args:

    Returns:
            a list of the percentage of similarity between our seq and the different genome
    """
    list_simili = [0] * nbr_colors
    index = 0
    for index in range(0,len(seq) - k):
        kmer = seq[index : index + k]
        for unitig in cdbg.keys():
            if kmer in unitig:
                for color in cdbg[unitig]:
                    list_simili[color]+=1
                break
    for i_list in range(nbr_colors):
        list_simili[i_list] = round(list_simili[i_list] / (len(seq) - k), 4)
    return list_simili

def find_colors(cdbg):
    colors = set()
    for c in cdbg.values():
        for val in c:
            colors.add(val)
    return len(colors)

def query_compute(file_name : str, cdbg : dict, k : int):
    nbr_colors = find_colors(cdbg)
    out = ""
    for record in SeqIO.parse(file_name,"fasta"):
        seq = str(record.seq)
        similarity = query_similarity(cdbg, seq, k, nbr_colors)
        out += f"{record.id} {similarity}\n"
    return(out)


# query_compute("Projet_alg/First_set/query.fa", {"a":"ACEHHE"})