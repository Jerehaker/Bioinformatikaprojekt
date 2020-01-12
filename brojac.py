f = open("test11tangle.gfa")
brojac = 0
brojac3 = 0
brojac2 = 0
while(True):
    line = f.readline()
    if not line:
        break

    if(line[0]=="S"):
        print(line.split("\t"))
        print(len(line.split("\t")[2])-1)
        brojac += len(line.split("\t")[2]) - 1
        brojac3 += 1
    if(line[0] == "L"):
        brojac2 += 1
print(brojac*2)
print(brojac2)
print(brojac3*4)
print(brojac*2 + brojac2)