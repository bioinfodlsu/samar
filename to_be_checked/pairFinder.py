import gc
import csv

p1 = open("paired_1_single.csv")
p2 = open("paired_2_single.csv")

#End Loop
Done = False

#Remove Header
p1.readline()
p2.readline()

# Initial lines
p1_line = p1.readline().rstrip()
p1_temp = p1_line.split(",")
p2_line = p2.readline().rstrip()
p2_temp = p2_line.split(",")

#new CSV header
final = [["Reference","Query_1","R_Index_1","Q_Index_1","R_Alen_1","Q_Alen_1","R_Len_1","Q_Len_1","Score_1","E_Value_1",
                      "Query_2","R_Index_2","Q_Index_2","R_Alen_2","Q_Alen_2","R_Len_2","Q_Len_2","Score_2","E_Value_2"]]
temp = [["","","","","","","","","","","","","","","","","","",""]]

#
count = 0

while not Done:
    #Compare read names
    if(p2_line == ""):
        p1_line = p1.readline().rstrip()
    elif(p1_line == ""):
        p2_line = p2.readline().rstrip()
    elif(p1_temp[1].split("/")[0] == p2_temp[1].split("/")[0] and p1_temp[0] == p2_temp[0]):
        for i in range(10):
            temp[0][i] = p1_temp[i]
        for i in range(1,10):
            temp[0][i+9] = p2_temp[i]
        final+= temp
        temp = [["","","","","","","","","","","","","","","","","","",""]]
        p2_line = p2.readline().rstrip()
        if(p2_line != ""):
            p2_temp = p2_line.split(",")
        p1_line = p1.readline().rstrip()
        if(p1_line != ""):
            p1_temp = p1_line.split(",")
    elif(int(p1_temp[1].split("_")[2].split("/")[0]) > int(p2_temp[1].split("_")[2].split("/")[0])):
        p2_line = p2.readline().rstrip()
        if(p2_line != ""):
            p2_temp = p2_line.split(",")
    elif(int(p1_temp[1].split("_")[2].split("/")[0]) < int(p2_temp[1].split("_")[2].split("/")[0])):
        p1_line = p1.readline().rstrip()
        if(p1_line != ""):
            p1_temp = p1_line.split(",")
    elif(int(p1_temp[1].split("_")[2].split("/")[0]) == int(p2_temp[1].split("_")[2].split("/")[0])):
        p2_line = p2.readline().rstrip()
        if(p2_line != ""):
            p2_temp = p2_line.split(",")
        p1_line = p1.readline().rstrip()
        if(p1_line != ""):
            p1_temp = p1_line.split(",")

    if(p1_line == "" and p2_line == ""):
        print("Done")
        Done = True

    gc.collect()
    count+=1
    
with open('paired_reads.csv', 'w', newline='') as file:
       writer = csv.writer(file)
       writer.writerows(final)

p1.close()
p2.close()
