
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
        if kmer in cdbg.keys():
            for color in cdbg[kmer]:
                list_simili[color]+=1
    for index in range(len(list_simili)):
        list_simili[index] = list_simili[index] / (len(seq) - k)
    return list_simili

def find_colors(cdbg):
    colors = set()
    for c in cdbg.values():
        for val in c:
            colors.add(val)
    return len(colors)

def query_compute(file_name : str, cdbg : dict):
    k = len(next(iter(cdbg)))
    nbr_colors = find_colors(cdbg)
    total_file = ""
    out = []
    with open(file_name,"r") as f:       #compute the diff seq into one str           
        for line in f:  
            line = line.strip()
            if ">" in line:
                line = "\n"+line+"\n"                          
            total_file += line
    for line in total_file.split(">"): #extract name and seq from the fasta str
        nom_seq = line.split("\n")[0]
        seq = line.split("\n")[1]
        similarity = query_similarity(cdbg, seq, k, nbr_colors)
        out.append((nom_seq,similarity))
    return(out)


# query_compute("Projet_alg/First_set/query.fa", {"a":"ACEHHE"})