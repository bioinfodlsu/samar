import csv

fp = open(snakemake.input[0])
Ref = True
final = [["Reference","Query","R_Index","Q_Index","R_Alen","Q_Alen","R_Len","Q_Len","Score","E_Value"]]
temp = [["","","","","","","","","",""]]
for i, line in enumerate(fp):
    if(line[0] == "a"):
    	tLine = line.split(" ")
    	temp[0][8] = tLine[1].split("=")[1]
    	temp[0][9] = tLine[3].split("=")[1]
    elif(line[0] == "s" and Ref):
        tLine = line.split()
        temp[0][0] = tLine[1]
        temp[0][2] = tLine[len(tLine)-5]
        temp[0][4] = tLine[len(tLine)-4]
        temp[0][6] = tLine[len(tLine)-2]
        Ref = False
    elif(line[0] == "s" and not Ref):
        tLine = line.split()
        temp[0][1] = tLine[1]
        temp[0][3] = tLine[len(tLine)-5]
        temp[0][5] = tLine[len(tLine)-4]
        temp[0][7] = tLine[len(tLine)-2]
        Ref = True
        final+= temp
        temp = [["","","","","","","","","",""]]

with open(snakemake.output[0], 'w', newline='') as file:
       writer = csv.writer(file)
       writer.writerows(final)

fp.close()
