#test
total_file = ""
total_file = ""
with open("query.fa","r") as f:                  
        for line in f:  
            line = line.strip()
            if ">" in line:
                line = "\n"+line+"\n"                          
            total_file += line
for line in total_file.split(">"):
    nom_seq = line.split("\n")[0]
    seq = line.split("\n")[1]
    print(nom_seq+"\n"+seq+'\n\n')

                   
                          