
def query_similarity(cdbg : dict, q : str , k : int):
    """
    Args:

    Returns:
            a list of the percentage of similarity between our seq and the different genome
    """
    list_simili = [0] * len(cdbg.values())
    for index in range(0,len(q) - k):
        kmer = q[index : index + k]
        if kmer in cdbg.keys():
            for color in cdbg[kmer]:
                list_simili[color]+=1
    for index in range(len(list_simili)):
        list_simili[index] = list_simili[index] / (len(q) - k)
    return list_simili

def query_compute(file_name : str, cdbg : dict):
    k = len(next(iter(cdbg)))
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
        similarity = query_similarity(cdbg,seq,k)
        out.append((nom_seq,similarity))
    print(out)


# query_compute("Projet_alg/First_set/query.fa", {"a":"ACEHHE"})