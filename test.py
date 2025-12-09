#test
# total_file = ""
# total_file = ""
# with open("query.fa","r") as f:                  
#         for line in f:  
#             line = line.strip()
#             if ">" in line:
#                 line = "\n"+line+"\n"                          
#             total_file += line
# for line in total_file.split(">"):
#     nom_seq = line.split("\n")[0]
#     seq = line.split("\n")[1]
#     print(nom_seq+"\n"+seq+'\n\n')

a = "AAAA"
dictio = {"GRTSUQDFUAAAAEHFSIFS":{1,2},"GRTSUQDFUAAAAAEHFSIFS":{1,4},"GRTSUQDFUAAEHFSIFS":{5,6}}
setiset = set()
for clé in dictio.keys():
    if a in clé:
        setiset.update(dictio[clé])
print(set(dictio.keys()))

                   
                          