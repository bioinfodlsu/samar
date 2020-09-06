import itertools
import csv
import argparse

def FileFilter(infile):
    for line in infile:
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

def main(input_alns, counts):

    final = {}
    lengths = {}
    median = 66
    iqr= 27
    upper = median + iqr
    lower = median - iqr

    def is_concordant(p1,p2):
        # 1:ref id, 9: query strand, 2: ref start pos 3: ref aln length
        if (p1[1] == p2[1] and p1[9] != p2[9]):
            #print("same chrom, diff strands")
            if p1[9] == '+' and int(p1[2]) < int(p2[2]):
                #print("directionality ok")
                frag_len = int(p2[2]) + int(p2[3]) - int(p1[2])
                #print(frag_len)
                return (lower <= frag_len <= upper)
            elif p2[9] == '+' and int(p2[2]) < int(p1[2]) :
                frag_len = int(p1[2]) + int(p1[3]) - int(p2[2])
                #print(frag_len)
                return (lower <= frag_len <= upper )
        #print("False")
        return False

    with open(input_alns) as infile:
        for key,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):
            block = list(group)
            #separate into 2 groups (pair 1 and pair 2)
            pair1 = [x for x in block if x[6] == key+"/1"]
            pair2 = [x for x in block if x[6] == key+"/2"]

            #count only concordant pairs for now
            concordant_pairs = []
            for p1 in pair1:
                for p2 in pair2:
                    #print(p1)
                    #print(p2)
                    if is_concordant(p1,p2) :
                        concordant_pairs.append(p1) #getting only one of the pairs

            num_of_corcondant_pairs = len(concordant_pairs)
            for aln in concordant_pairs:
                ref_id = p1[1]
                ref_len = p1[5]
                final[ref_id] = final.get(ref_id,0) + 1.0/num_of_corcondant_pairs
                if not ref_id in lengths: lengths[ref_id] = ref_len
            #print("====")
            # #Case 1: pair 1 and 2 have exactly 1 alignment
            # #if((len(pair1) == 1 and len(pair2) == 1) or (len(pair1) == 1 and len(pair2) == 0) or (len(pair1) == 0 and len(pair2) == 1)):
            # if(len(pair1) == 1 and len(pair2) == 1):
            #     for alg in block:
            #         final[alg[1]] = final.get(alg[1],0) + 1
            # #Case 2: pair 1 or pair 2 has more than 1 alignment and both are not len 0
            # elif(len(pair1) >= 1 and len(pair2) >= 1):
            #     temp = []
            #     for p1 in pair1:
            #         for p2 in pair2:
            #             if(p1[1] == p2[1] and lower <= GetAlignment(p1,p2) <= upper):
            #                 temp.append(p1[1])
            #
            #     for ref in temp:
            #         final[ref] = final.get(ref,0) + (1/len(temp))


    with open(counts,"w") as tab_file:
        writer = csv.writer(tab_file, delimiter='\t')

        #adding Header
        writer.writerow(["Ref_Seq_ID","Seq_length","Num_of_Alg"])

        for key, value in final.items():
            writer.writerow([key,lengths[key],value])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',help="alignments in tab format")
    parser.add_argument('output_file',help="output file containing counts")
    args = parser.parse_args()
    main(args.input_file,args.output_file)
