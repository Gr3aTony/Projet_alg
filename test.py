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
dictio.update({"aezfz":{1,3,4}})
print(dictio)
lis = ["zlfe,ze",2,3]
sete = {1,2,3,6}
sete.difference_update(lis)
print(sete)
sete.add(lis[0])

az ="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
print(az[:-(len(az)-1)])
a = set()                 
# sete.discard(1)
a = a.union(sete)
a.discard(1)
print(a)