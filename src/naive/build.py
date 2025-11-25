def dbg_builder(fasta_file,color, k,dbg): # fasta_file : str || color : int || k : int || dbg : dict
    """
    Creer le dictionnaire contenant le graphe de De Brujns coloré
    fasta_file = pwd du fichier fasta,color indice du fichier dans sa liste, k taille des kmers,dbg dictionaire du graphe
    renvoi le dictionaire mis à jours
    """
    with open(fasta_file,"r") as f:                  # Ouverture du fichier FASTA
        for line in f:                               
            line = line.strip()                 # Suppression des caractères spéciaux et espaces pour n'avoir (en théorie) que ATCG.
            if not line.startswith('>'):        # On ignore les lignes de description FASTA
                for i in range(0,len(line)-k):
                    kmer = line[i:i+k]          # A chaque ittération on avance de 1 dans la séquence en crée un kmer temporaire composé du caractère à l'indice actuel + les k caractères suivant
                    # Ajout du k-mer et de sa couleur
                    if kmer not in dbg.keys(): #Si le kmer n'existe pas encore il est crée et la couleur actuelle y est ajoutée
                        dbg[kmer] = [color]
                    elif color not in dbg[kmer]: # Si la couleur n'est pas attribué à ce kmer dans le dictionnaire on l'y attribue
                        dbg[kmer].append(color)
    return dbg



def loop(file_list,k):
    dbg = {}
    i = 0 
    with open(file_list,"r") as fl:                  # Ouverture du fichier contenant les chemins d'accès aux fichiers génomiques
        for file in fl:
            file =file.strip()                        # Suppression des caractères spéciaux et espaces pour n'avoir que les chemins des fichiers.
            dbg = dbg_builder(file,i,k,dbg)
            i+=1
    return dbg

