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


def right_unitig(kmer,dbg):
    unitig = kmer
    current_kmer = kmer
    while True:
        voisins = right_neighbors(current_kmer,dbg)
        if len(voisins) == 0:
            print("Stop because dead-end")
            return unitig
        if len(voisins) >1:
            print("Stop because multiple out edges")
            return unitig
        if voisins[0] == kmer:
            print("Stop because loop")
            return unitig
        if len(left_neighbors(voisins[0],dbg))>1:
            print("Stop because multiple in edges")
            return unitig
        # else can extend
        unitig += voisins[0][-1]
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

    return dbg,bag_of_kmer
