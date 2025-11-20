def dbg_builder(fasta_file,color, k,dbg):
    with open(fasta_file,"r") as f:
        for line in f:
            line = line.strip()
            print(line)
            if not line.startswith('>'):
                for i in range(0,len(line)-k):
                    kmer = line[i:i+k]
                    if kmer not in dbg.keys():
                        dbg[kmer] = [color]
                    else:
                        dbg[kmer].append(color)
    return dbg

