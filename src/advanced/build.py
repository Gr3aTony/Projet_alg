from Bio import SeqIO
def dbg_builder(fasta_file : str, color : int, k : int, dbg : dict): # fasta_file : str || color : int || k : int || dbg : dict
    """
    Updates De Brujns colored graph of a single fasta file
    Args:
        fasta_file : pwd of a fasta file ,
        color : place of the file in the original list
        k : size of kmers
        dbg : De brujn graph's 
    Returns:
        dict : an updated De Brujns colored graph
    
    """
    for record in SeqIO.parse(fasta_file,"fasta"):
        seq = str(record.seq)   # On ignore les lignes de description FASTA
        for i in range(0, len(seq) - k):
            kmer = seq[i : i + k]          # A chaque ittération on avance de 1 dans la séquence en crée un kmer temporaire composé du caractère à l'indice actuel + les k caractères suivant
            # Ajout du k-mer et de sa couleur
            if kmer not in dbg.keys(): #Si le kmer n'existe pas encore il est crée et la couleur actuelle y est ajoutée
                dbg[kmer] = [color]
            elif color not in dbg[kmer]: # Si la couleur n'est pas attribué à ce kmer dans le dictionnaire on l'y attribue
                dbg[kmer].append(color)
    return dbg

def right_neighbors(kmer, dbg):
    res = []
    for letter in ["A", "C", "G", "T"]:
        voisin = kmer[1:] + letter
        if voisin in dbg:
            res.append(voisin)
    return res

def left_neighbors(kmer, dbg):
    res = []
    for letter in ["A", "C", "G", "T"]:
        voisin = letter + kmer[:-1]
        if voisin in dbg:
            res.append(voisin)
    return res

def right_unitig(kmer,list_kmer):
    unitig = kmer
    used = {kmer}
    current = kmer

    while True:
        voisins = right_neighbors(current, list_kmer)
        if len(voisins) != 1:              
            return unitig, used

        nxt = voisins[0]
        unitig += nxt[-1]               
        used.add(nxt)

        if nxt == kmer or nxt == current:
            return unitig, used
        if len(left_neighbors(nxt, list_kmer)) != 1 or len(right_neighbors(nxt, list_kmer)) != 1:
            return unitig, used

        current = nxt

def left_unitig(kmer,list_kmer):
    unitig = kmer
    used = {kmer}
    current = kmer

    while True:
        voisin = left_neighbors(current, list_kmer)
        if len(voisin) != 1:              
            return unitig, used

        nxt = voisin[0]
        unitig = nxt[0] + unitig        
        used.add(nxt)

        if nxt == kmer or nxt == current:#en cas de boucle ou de genome circulaire
            return unitig, used
        if len(left_neighbors(nxt, list_kmer)) != 1 or len(right_neighbors(nxt, list_kmer)) != 1:
            return unitig, used

        current = nxt
        


def loop(file_list : str, k : int):
    """
    Args :
            file_list: list of pwd of multiple fasta files
            k : size of kmers
    """
    dbg = {}
    i = 0 
    with open(file_list,"r") as fl:   # Ouverture du fichier contenant les chemins d'accès aux fichiers génomiques
        for file in fl:
            file = file.strip()     # Suppression des caractères spéciaux et espaces pour n'avoir que les chemins des fichiers.
            dbg = dbg_builder(file, i, k, dbg)
            i += 1
    bag_of_kmer = set(dbg.keys())
    list_kmer = set(dbg.keys())
    final_dict = {}

    while list_kmer:
        current = next(iter(list_kmer))
        colors_unitig = dbg[current]#recup couleur de l'unitig que l'on va creer

        left_uni, left_used = left_unitig(current, bag_of_kmer)
        right_uni, right_used = right_unitig(current, bag_of_kmer)

        unitig = left_uni[:-k] + right_uni
        used = left_used.union(right_used)
        final_dict[unitig] = sorted(colors_unitig)

        list_kmer.difference_update(used)

    return final_dict, bag_of_kmer
            



# dico,kmer = loop("G.txt",31)
# i=0
# for x in dico.keys():
#     for j in dico.keys():
#         if x in j and x!=j:
#             print(x,j.find(x),len(j)-len(x))
#             i+=1
# print(len(dico.keys()))
