import itertools
import csv

def FileFilter(file):
    for line in file:
        if not line.startswith("#"):
            yield line.split()

def GetAlignment(p1,p2):
    temp = 0
    
    if (int(p1[2]) > int(p2[2])):
        temp = int(p1[2]) - int(p2[2]) + int(p1[3])
    elif(int(p1[2]) < int(p2[2])):
        temp = int(p2[2]) - int(p1[2]) + int(p2[3])
    else:
        temp = int(p1[3])

    return temp

final = {}
median = 66
iqr= 27
upper = median + iqr
lower = median - iqr

with open("interleaved_alg.tab") as file:
    for key,group in itertools.groupby(FileFilter(file), lambda x : x[6].split("/")[0]):
        block = list(group)
        #separate into 2 groups (pair 1 and pair 2)
        pair1 = [x for x in block if x[6] == key+"/1"]
        pair2 = [x for x in block if x[6] == key+"/2"]
        #Case 1: pair 1 and 2 have exactly 1 alignment
        #if((len(pair1) == 1 and len(pair2) == 1) or (len(pair1) == 1 and len(pair2) == 0) or (len(pair1) == 0 and len(pair2) == 1)):
        if(len(pair1) == 1 and len(pair2) == 1):
            for alg in block:
                final[alg[1]] = final.get(alg[1],0) + 1
        #Case 2: pair 1 or pair 2 has more than 1 alignment and both are not len 0
        elif(len(pair1) >= 1 and len(pair2) >= 1):
            temp = []
            for p1 in pair1:
                for p2 in pair2:
                    if(p1[1] == p2[1] and lower <= GetAlignment(p1,p2) <= upper):
                        temp.append(p1[1])

            for ref in temp:
                final[ref] = final.get(ref,0) + (1/len(temp))

        #Currently Ignores Case if only 1 pair has more than 1 alignment 

with open("interleaved_alg_count.tab","w") as tab_file:
    writer = csv.writer(tab_file, delimiter='\t')
    
    #adding Header
    writer.writerow(["Ref_Seq","Num_of_Alg"]) 

    for key, value in final.items():
        writer.writerow([key,value])
    
