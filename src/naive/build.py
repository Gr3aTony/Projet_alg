def dbg_builder(fasta_file,color, k,dbg):
    """
    Creer le dictionnaire contenant le graphe de De Brujns coloré
    fasta_file = pwd du fichier fasta,color indice du fichier dans sa liste, k taille des kmers,dbg dictionaire du graphe
    renvoi le dictionaire mis à jours
    """
    with open(fasta_file,"r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                for i in range(0,len(line)-k):
                    kmer = line[i:i+k]
                    if kmer not in dbg.keys():
                        dbg[kmer] = [color]
                    elif color not in dbg[kmer]:
                        dbg[kmer].append(color)
    return dbg



def loop(file_list,k):
    dbg = {}
    i = 0 
    with open(file_list,"r") as fl:
        for file in fl:
            file =file.strip()
            dbg = dbg_builder(file,i,k,dbg)
            i+=1
    return dbg

print(loop("/Users/2ndlife/Documents/Rennes M2/ALG/Projet_alg/G.txt",17))