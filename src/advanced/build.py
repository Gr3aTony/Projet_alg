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

# def right_unitig(kmer,list_kmer):
#     unitig = kmer
#     current_kmer = kmer
#     list_neighbors = [kmer]
#     while True:
#         left_conect = left_neighbors(kmer,list_kmer)
#         voisins = right_neighbors(current_kmer,list_kmer)
#         if len(left_conect) == 1:
#             return 0
#         if len(voisins) == 0:
#             return unitig,list_neighbors
#         if len(voisins) >1:
#             return unitig,list_neighbors
#         if voisins[0] == kmer:
#             return unitig,list_neighbors
#         if len(left_neighbors(voisins[0],list_kmer))>1:
#             return unitig,list_neighbors
#         unitig += voisins[0][-1]
#         current_kmer = voisins[0]
#         list_neighbors.append(voisins[0])

# def left_unitig(kmer,list_kmer):
#     unitig = kmer
#     current_kmer = kmer
#     list_neighbors = [kmer]
#     while True:
#         right_connect = right_neighbors(kmer,list_kmer)
#         voisins = left_neighbors(current_kmer,list_kmer)
#         if len(right_connect)>1:
#             return 0
#         if len(voisins) == 0:
#             return unitig,list_neighbors
#         if len(voisins) >1:
#             return unitig,list_neighbors
#         if len(right_neighbors(voisins[0],list_kmer))>1:
#             return unitig,list_neighbors

#         unitig = voisins[0][0]+unitig
#         current_kmer = voisins[0]
#         list_neighbors.append(voisins[0])

def nbr_voisin_gauche(kmer, kmer_set):
    return len(left_neighbors(kmer, kmer_set))
def nbr_voisin_droit(kmer, kmer_set):
    return len(right_neighbors(kmer, kmer_set))

def find_start(kmer, kmer_set):
    while True:
        if not (nbr_voisin_gauche(kmer, kmer_set) == 1 ):
            return kmer
        left = left_neighbors(kmer, kmer_set)[0]
        if nbr_voisin_droit(left, kmer_set) != 1:
            return kmer
        kmer = left

def build_unitig(start, kmer_set):
    unitig = start
    used = {start}
    current = start

    while nbr_voisin_droit(current, kmer_set) == 1:
        nxt = right_neighbors(current, kmer_set)[0]
        if nxt in used:          # cycle
            break
        used.add(nxt)
        unitig += nxt[-1]
        current = nxt
        # stop after adding an endpoint (endpoint = not internal)
        if not (nbr_voisin_gauche(current, kmer_set) == 1 and nbr_voisin_droit(current, kmer_set) == 1):
            break

    return unitig, used


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
        colors_unit = dbg[current]

        start_uni = find_start(current,list_kmer)
        indeg = len(left_neighbors(start_uni, list_kmer))
        outdeg = len(right_neighbors(start_uni, list_kmer))
        assert (indeg==1 and outdeg==1), (start_uni, indeg, outdeg)
        unitig,used = build_unitig(start_uni,list_kmer)
        final_dict[unitig] =colors_unit
        list_kmer.difference_update(used)
            
    return final_dict,bag_of_kmer


test = loop("G.txt",31)
print(len(test[0]))
