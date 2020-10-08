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

def main(input_alns, out_counts, multi_mapping, frag_len_mean, frag_len_std):

    final = {}
    lengths = {}


    if multi_mapping:
        '''
        1. Because we are only aligning to the (translated) CDS, there could be many cases with only one end aligning. Treat them as if they are one fragment.
        2. If we do see a concodarnt pair, it's most likely that any other unpaired alignments are of low quality and safe to ignore.
        '''

        lower = frag_len_mean - 3* frag_len_std
        upper = frag_len_mean + 3* frag_len_std
        
        def update_from_single_end(alns):
            for aln in alns:
                ref_id = aln[1]
                ref_len = aln[5]
                mismap = float(aln[-1].split("=")[-1])
                final[ref_id] = final.get(ref_id,0) + (1-mismap)
                if not ref_id in lengths : lengths[ref_id] = ref_len

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
                    return (lower <= frag_len <= upper)
            #print("False")
            return False


        with open(input_alns) as infile:
            for key,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):

                block = list(group)

                alns1 = [x for x in block if x[6] == key+"/1"] #alignments of read 1
                alns2 = [x for x in block if x[6] == key+"/2"] #aligments of read 2

                if len(alns2) == 0:
                    update_from_single_end(alns1)
                elif len(alns1) == 0:
                    update_from_single_end(alns2)
                else:
                    #get pairs. ignore singletons. distribute counts using mismap of one of the reads.
                    for p1 in alns1:
                        for p2 in alns2:
                            if is_concordant(p1,p2):
                                ref_id = p1[1]
                                ref_len = p1[5]
                                mismap = float(p1[-1].split("=")[-1])
                                final[ref_id] = final.get(ref_id,0) + (1-mismap)
                                if not ref_id in lengths: lengths[ref_id] = ref_len

    else:
        #the assumption here is that there is at most 1 alignment per read,
        #i.e. last-pair-probs was run with -m less than 0.5

        with open(input_alns) as infile:
            for pair_id,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):

                block = list(group)

                alns1 = [x for x in block if x[6] == pair_id+"/1"]
                alns2 = [x for x in block if x[6] == pair_id+"/2"]

                try:
                    assert(len(alns1) <= 1)
                    assert(len(alns2) <= 1)
                except AssertionError:
                    print("Too many alignments for read pair {}".format(pair_id))
                    continue

                if len(alns1) ==  1 and len(alns2) == 1 :
                    if alns1[0][9] == alns2[0][9]:
                        final[alns1[0][1]] = final.get(alns1[0][1],0) + 1
                        if not alns1[0][1] in lengths : lengths[alns1[0][1]] = alns1[0][5]

                elif len(alns1) == 1:
                    final[alns1[0][1]] = final.get(alns1[0][1],0) + 1
                    if not alns1[0][1] in lengths : lengths[alns1[0][1]] = alns1[0][5]

                elif len(alns2) == 1:
                    final[alns2[0][1]] = final.get(alns2[0][1],0) + 1
                    if not alns2[0][1] in lengths : lengths[alns2[0][1]] = alns2[0][5]



    with open(out_counts,"w") as tab_file:
        writer = csv.writer(tab_file, delimiter='\t')

        #adding Header
        writer.writerow(["Ref_Seq_ID","Seq_length","Num_of_Alg"])

        for key, value in final.items():
            writer.writerow([key,lengths[key],value])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',help="alignments in tab format")
    parser.add_argument('output_file',help="output file containing counts")
    parser.add_argument('--multi_mapping',action='store_true')
    parser.add_argument('--frag_len_mean',type=float,help="mean of fragment length, when translated. Use the same value as last-pair-probs")
    parser.add_argument('--frag_len_std',type=float,help="standard deviation of fragment length, when translated. Use the same value as last-pair-probs")

    args = parser.parse_args()
    main(args.input_file,args.output_file,args.multi_mapping, args.frag_len_mean, args.frag_len_std)
