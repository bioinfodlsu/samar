#%%
import itertools
import csv
import argparse
import copy

def FileFilter(infile):
    for line in infile:
        if not line.startswith("#"):
            yield line.split()

#%%
def main(input_alns, out_counts, frag_len_mean, frag_len_std):
    '''
    Does 2 passes over the alignments. 
    In the first pass, we only consider reads with unique alignments. Counts are recorded in the dict unique.
    In the second pass, we consider the remaining reads. Counts are distributed based on the proportion of uniquely mapped reads. 
    For any read pair, if both reads have 
    '''
    
    lower = frag_len_mean - 3* frag_len_std
    upper = frag_len_mean + 3* frag_len_std
    #final = {}
    unique = {} #key = peptide ID, value = uniquely mapped fragment counts
    lengths = {}
 
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
 
    #First pass
    with open(input_alns) as infile:
        for key,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):

            block = list(group)

            alns1 = [x for x in block if x[6] == key+"/1"] #alignments of read 1
            alns2 = [x for x in block if x[6] == key+"/2"] #alignments of read 2
            
            #initialize unique[ref_id] to 0. need to do this outside the unique checking logic so that we are ready for second pass
            for aln in alns1:
                unique[aln[1]] = unique.get(aln[1],0)
                if not aln[1] in lengths: lengths[aln[1]] = aln[5]
            
            for aln in alns2:
                unique[aln[1]] = unique.get(aln[1],0)
                if not aln[1] in lengths: lengths[aln[1]] = aln[5]
            
            #if its a unique pair
            if len(alns1) ==  1 and len(alns2) == 1 : 
                if is_concordant(alns1[0], alns2[0]):
                   unique[alns1[0][1]] += 1
                   
            #if only read1 has an alignment and its unique
            elif len(alns1) == 1 and len(alns2) == 0:
                unique[alns1[0][1]] += 1

            #if only read2 has an alignment and its unique
            elif len(alns2) == 1 and len(alns1) == 0:
                unique[alns2[0][1]] += 1
                

    #second pass
    final = copy.deepcopy(unique)
        
   
    def update_from_single_end(alns):
        ref_ids = [aln[1] for aln in alns]
        unique_counts = [unique[ref_id] for ref_id in ref_ids]
        if sum(unique_counts) != 0 :
            props = [c/sum(unique_counts) for c in unique_counts]
            for i,ref_id in enumerate(ref_ids):
                final[ref_id] += props[i] 
        else:
            denom = len(ref_ids)
            for ref_id in ref_ids:
                final[ref_id] += 1/denom
        
    with open(input_alns) as infile:
        for key,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):

            block = list(group)

            alns1 = [x for x in block if x[6] == key+"/1"] #alignments of read 1
            alns2 = [x for x in block if x[6] == key+"/2"] #alignments of read 2
            
            if len(alns1) > 1 and len(alns2) > 1: 
                ref_ids = [] #counts need to be updated for these
                for p1 in alns1:
                    for p2 in alns2:
                        if is_concordant(p1,p2):
                            ref_ids.append(p1[1])
                unique_counts = [unique[ref_id] for ref_id in ref_ids]
                if sum(unique_counts) != 0 :
                    props = [c/sum(unique_counts) for c in unique_counts]
                    for i,ref_id in enumerate(ref_ids):
                        final[ref_id] += props[i] 
                else: #distribute evenly
                    denom = len(ref_ids)
                    for ref_id in ref_ids:
                        final[ref_id] += 1/denom
                    
            elif len(alns1) > 1 and len(alns2) == 0: 
                update_from_single_end(alns1)
                
            elif len(alns2) > 1 and len(alns1) ==0:
                update_from_single_end(alns2)
            
            else:
                continue
            
    #unique_c = [unique[k] for k in unique.keys()]
    #final_c = [final[k] for k in unique.keys()]
    #return (unique_c,final_c)

    with open(out_counts,"w") as tab_file:
        writer = csv.writer(tab_file, delimiter='\t')

        #adding Header
        writer.writerow(["Ref_Seq_ID","Seq_length","Num_of_Alg"])

        for key, value in final.items():
            writer.writerow([key,lengths[key],value])
    
    # for k,v in final.items():
    #     print('{0} :{1}'.format(k,v))


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
# u,f = main("/home/anish/Desktop/last_multimap/1mil.tab","out",82,8)
# import matplotlib.pyplot as plt
# import numpy as np

# u_1 = [ item for item in u if item < 10 ]
# f_1 = [ item for item in f if item < 10 ]

# plt.hist([u_1,f_1],bins=10, label=['u','f'])
