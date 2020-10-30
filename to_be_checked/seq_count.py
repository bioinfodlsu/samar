#%%
import itertools
import csv
import argparse
import copy

class Tab_alns:
    def __init__(self,tab_line):
        self.ref = tab_line[1]
        self.ref_start = int(tab_line[2])
        self.ref_aln_width = int(tab_line[3])
        self.ref_length = int(tab_line[5])
        self.query = tab_line[6]
        self.query_strand = tab_line[9]


def FileFilter(infile):
    for line in infile:
        if not line.startswith("#"):
            yield line.split()

def is_concordant(p1,p2): #given two alignments, returns (True,start,end) if concordant and (False,0,0) otherwise
    # 1:ref id, 9: query strand, 2: ref start pos 3: ref aln length
    if (p1.ref == p2.ref and p1.query_strand != p2.query_strand): #same ref seq and different strands

        if p1.query_strand == '+' and p1.ref_start < p2.ref_start:
            if (lower <= p2.ref_start + p2.ref_aln_width - p1.ref_start <= upper):
                return (True, p1.ref_start, p2.ref_start + p2.ref_aln_width-1)

        elif p2.query_strand == '+' and p2.ref_start < p1.ref_start:
            if (lower <= p1.ref_start + p1.ref_aln_width - p2.ref_start <= upper):
                return (True, p2.ref_start, p1.ref_start + p1.ref_aln_width-1 )

    return(False,0,0)

def unique_pass(input_alns,unique):

    def update_unique_single_end(aln):
        unique[aln.ref][-1] += 1 #update count
        for i in range(aln.ref_start,aln.ref_start + aln.ref_aln_width):
            unique[aln.ref][i] +=1


    with open(input_alns) as infile:
        for key,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):

            block = list(group)
            alns1 = [Tab_alns(x) for x in block if x[6] == key+"/1"] #alignments of read 1
            alns2 = [Tab_alns(x) for x in block if x[6] == key+"/2"] #alignments of read 2

            #initialize unique[ref_id] to a 0-vector of length of ref_id sequence. need to do this outside the unique checking logic so that we are ready for second pass
            for aln in alns1+alns2:
                unique[aln.ref] = unique.get(aln.ref,[0]*(aln.ref_length+1)) #final element is the count of fragments


            #if its a unique pair
            if len(alns1) ==  1 and len(alns2) == 1 :
                conc,frag_start,frag_end = is_concordant(alns1[0], alns2[0])
                if conc:
                   unique[alns1[0].ref][-1] += 1 #update count
                   for i in range(frag_start,frag_end+1): #assuming gapless aln
                       unique[alns1[0].ref][i] +=1

            #if only read1 has an alignment and it's unique
            elif len(alns1) == 1 and len(alns2) == 0:
                update_unique_single_end(alns1[0])


            #if only read2 has an alignment and it's unique
            elif len(alns2) == 1 and len(alns1) == 0:
                update_unique_single_end(alns2[0])

    #normalize unique counts by length (of non-zero counts)
    for k,v in unique.items():
        nonZeros = len(v)-1-v[:-1].count(0)
        if v[-1] != 0 :
            v.append(v[-1]/nonZeros)
        else:
            v.append(0)

def rescue_pass(input_alns,unique, final):

    def update_rescue_single_end(alns):
        ref_ids = [aln.ref for aln in alns]
        unique_norm_counts = [unique[ref_id][-1] for ref_id in ref_ids]
        denom = sum(unique_norm_counts)
        if denom != 0 :
            props = [c/denom for c in unique_norm_counts]
            for i,ref_id in enumerate(ref_ids):
                final[ref_id] += props[i]
        else:
            denom = len(ref_ids)
            for ref_id in ref_ids:
                final[ref_id] += 1/denom

    with open(input_alns) as infile:
        for key,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):

            block = list(group)

            alns1 = [Tab_alns(x) for x in block if x[6] == key+"/1"] #alignments of read 1
            alns2 = [Tab_alns(x) for x in block if x[6] == key+"/2"] #alignments of read 2

            if len(alns1) > 1 and len(alns2) > 1:
                ref_ids = [] #counts need to be updated for these
                for p1 in alns1:
                    for p2 in alns2:
                        if is_concordant(p1,p2):
                            ref_ids.append(p1.ref)

                unique_norm_counts = [unique[ref_id][-1] for ref_id in ref_ids]
                denom = sum(unique_norm_counts)
                if denom != 0 :
                    props = [c/denom for c in unique_norm_counts]
                    for i,ref_id in enumerate(ref_ids):
                        final[ref_id] += props[i]
                else: #distribute evenly
                    denom = len(ref_ids)
                    for ref_id in ref_ids:
                        final[ref_id] += 1/denom

            elif len(alns1) > 1 and len(alns2) == 0:
                update_rescue_single_end(alns1)

            elif len(alns2) > 1 and len(alns1) == 0:
                update_rescue_single_end(alns2)

            else:
                continue


def main(input_alns, out_counts, frag_len_mean, frag_len_std):
    '''
    Does 2 passes over the alignments.
    In the first pass, we only consider reads with unique alignments. Counts are recorded in the dict unique.
    In the second pass, we consider the remaining reads. Counts are distributed based on the proportion of uniquely mapped reads.
    For any read pair, if both reads have
    '''

    global lower, upper
    lower = frag_len_mean - 3* frag_len_std
    upper = frag_len_mean + 3* frag_len_std

    unique = {} #key = peptide ID, value = list of length size of peptide with number of fragments covering each position

    #first pass
    unique_pass(input_alns,unique)

    #second pass
    global final
    final={}
    for k,v in unique.items():
        final[k] = v[-2] #current unique counts
    #final = copy.deepcopy(unique)
    rescue_pass(input_alns,unique, final)

    #output
    with open(out_counts,"w") as tab_file:
        writer = csv.writer(tab_file, delimiter='\t')

        #adding Header
        writer.writerow(["Ref_Seq_ID","Seq_length","Num_of_Alg"])

        for key, value in final.items():
            writer.writerow([key,len(unique[key])-2,value])
    return (unique,final)

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',help="alignments in tab format")
    parser.add_argument('output_file',help="output file containing counts")
    #parser.add_argument('--multi_mapping',action='store_true')
    parser.add_argument('--frag_len_mean',type=float,help="mean of fragment length, when translated. Use the same value as last-pair-probs")
    parser.add_argument('--frag_len_std',type=float,help="standard deviation of fragment length, when translated. Use the same value as last-pair-probs")

    args = parser.parse_args()
    main(args.input_file,args.output_file, args.frag_len_mean, args.frag_len_std)

#%%
# unique,final=main("/home/anish/Desktop/last_multimap/1mil.tab","out",82,8)
#
# uf={}
# for k,v in final.items():
#      uf[k] = (unique[k][-2],v)
#
#
# import matplotlib.pyplot as plt
# import numpy as np
#
# u_1 = [ item for item in u if item < 10 ]
# f_1 = [ item for item in f if item < 10 ]
#
# plt.hist([u_1,f_1],bins=10, label=['u','f'])
