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
    current_kmer = kmer
    used_kmer = set()
    while True:
        voisins = right_neighbors(current_kmer,list_kmer)
        if len(voisins) == 0:
            used_kmer.add(current_kmer)
            return unitig,used_kmer
        if len(voisins) >1:
            return unitig,used_kmer
        if voisins[0] == kmer:
            used_kmer.add(current_kmer)
            return unitig,used_kmer
        if len(left_neighbors(voisins[0],list_kmer))>1:
            used_kmer.add(current_kmer)
            return unitig,used_kmer
        unitig += voisins[0][-1]
        used_kmer.add(current_kmer)
        current_kmer = voisins[0]
        

def left_unitig(kmer,list_kmer):
    unitig = kmer
    current_kmer = kmer
    used_kmer = set()
    while True:
        voisins = left_neighbors(current_kmer,list_kmer)
        if len(voisins) == 0:
            used_kmer.add(current_kmer)
            return unitig,used_kmer
        if len(voisins) >1:
            used_kmer.add(current_kmer)
            return unitig,used_kmer
        if voisins[0] == kmer:
            used_kmer.add(current_kmer)
            return unitig,used_kmer
        if len(right_neighbors(voisins[0],list_kmer))>1:
            used_kmer.add(current_kmer)
            return unitig,used_kmer
        unitig = voisins[0][0]+unitig
        used_kmer.add(current_kmer)
        current_kmer = voisins[0]
        


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
    seen = set()
    while list_kmer!=seen and len(list_kmer)!=0:
        # print(len(list_kmer))
        candidat = list_kmer.difference(seen)#tt kmer pas encore analysé
        current = next(iter(candidat))
        colors_unit = dbg[current]#recup couleur de l'unitig que l'on va creer
        left_uni,left_used = left_unitig(current,list_kmer)
        right_uni,right_used = right_unitig(current,list_kmer)
        if left_uni == right_uni:#cas ou on a un kmer unique avec 2 voisins a gauche et a droite
            used = left_used.union(right_used)
            list_kmer.difference_update(used)
            final_dict[left_uni] = colors_unit
        else:
            unitig = left_uni[:-k] + right_uni
            final_dict[unitig] =colors_unit
        if len(left_used) + len(right_used) <= 1:#current à un voisinage multiple a gauche ou a droite
            seen.add(current)
        elif len(left_used) == 0 or len(right_used) == 0 :# current a au moins un voisinage multiple direct
            used = left_used.union(right_used)
            # used.discard(current)
            list_kmer.difference_update(used)
        else:
            used = left_used.union(right_used)
            list_kmer.difference_update(used)
            
    return final_dict,bag_of_kmer


# dico,kmer = loop("G.txt",31)
# # i=0
# # for x in kmer:
# #     for j in dico.keys():
# #         if x in j:
# #             # print(x,j.find(x),len(j)-len(x))
# #             i+=1
# #             break
# print(len(dico.keys()))
