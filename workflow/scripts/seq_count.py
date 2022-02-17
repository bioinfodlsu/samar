#%%
import itertools
import csv
import argparse
from Bio import SeqIO

class Tab_alns:
    def __init__(self,tab_line):
        self.ref = tab_line[1]
        self.ref_start = int(tab_line[2])
        self.ref_aln_width = int(tab_line[3])
        self.ref_length = int(tab_line[5])
        self.query = tab_line[6]
        self.query_strand = tab_line[9]

class Counts:
    def  __init__(self,ref_length): 
        self.length = ref_length
        self.count_array = [0]*ref_length
        self.unique_count = 0.0
        self.unique_count_norm = 0.0
        self.final_count = 0.0
        self.final_count_norm = 0.0
        self.tpm = 0.0
    
    def __str__(self):
        return 'length:{self.length}    unique_count:{self.unique_count}    unique_count_norm:{self.unique_count_norm}   final_count:{self.final_count}    final_count_norm:{self.final_count_norm}    tpm:{self.tpm}'.format(self=self)

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

def unique_pass(input_alns,counts_dict):

    def update_unique_single_end(aln):
        counts_dict[aln.ref].unique_count += 1 #update count
        for i in range(aln.ref_start,aln.ref_start + aln.ref_aln_width):
            counts_dict[aln.ref].count_array[i] +=1


    with open(input_alns) as infile:
        for key,group in itertools.groupby(FileFilter(infile), lambda x : x[6].rsplit("/",1)[0]):

            block = list(group)
            alns1 = [Tab_alns(x) for x in block if x[6] == key+"/1"] #alignments of read 1
            alns2 = [Tab_alns(x) for x in block if x[6] == key+"/2"] #alignments of read 2

            
            #if it's a unique pair
            if len(alns1) ==  1 and len(alns2) == 1 :
                conc,frag_start,frag_end = is_concordant(alns1[0], alns2[0])
                if conc:
                   counts_dict[alns1[0].ref].unique_count += 1
                   for i in range(frag_start,frag_end+1): #assuming gapless aln
                       #unique[alns1[0].ref][i] +=1
                       counts_dict[alns1[0].ref].count_array[i] +=1


            #if only read1 has an alignment and it's unique
            elif len(alns1) == 1 and len(alns2) == 0:
                update_unique_single_end(alns1[0])


            #if only read2 has an alignment and it's unique
            elif len(alns2) == 1 and len(alns1) == 0:
                update_unique_single_end(alns2[0])

    #done with first pass. 1. set final_count to be unique_count and 2.normalize unique counts by length (of non-zero counts)
    for k,v in counts_dict.items():
	        
        v.final_count = v.unique_count
        
        #nonZeros = len(v)-1-v[:-1].count(0)
        nonZeros = len(v.count_array) - v.count_array.count(0)
        if v.unique_count != 0 :
            v.unique_count_norm = v.unique_count/nonZeros
        else:
            v.unique_count_norm = 0

def rescue_pass(input_alns,counts_dict):

    def update_rescue_single_end(alns):
        ref_ids = [aln.ref for aln in alns]
        unique_norm_counts = [counts_dict[ref_id].unique_count_norm for ref_id in ref_ids]
        denom = sum(unique_norm_counts)
        if denom != 0 :
            props = [c/denom for c in unique_norm_counts]
            for i,ref_id in enumerate(ref_ids):
                counts_dict[ref_id].final_count += props[i]
        else:
            denom = len(ref_ids)
            for ref_id in ref_ids:
                counts_dict[ref_id].final_count += 1/denom

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

                unique_norm_counts = [counts_dict[ref_id].unique_count_norm for ref_id in ref_ids]
                denom = sum(unique_norm_counts)
                if denom != 0 :
                    props = [c/denom for c in unique_norm_counts]
                    for i,ref_id in enumerate(ref_ids):
                        counts_dict[ref_id].final_count += props[i]
                else: #distribute evenly
                    denom = len(ref_ids)
                    for ref_id in ref_ids:
                        counts_dict[ref_id].final_count += 1/denom

            elif len(alns1) > 1 and len(alns2) == 0:
                update_rescue_single_end(alns1)

            elif len(alns2) > 1 and len(alns1) == 0:
                update_rescue_single_end(alns2)

            else:
                continue


def main(input_alns, out_counts, frag_len_mean, frag_len_std,reference):
    '''
    Does 2 passes over the alignments.
    In the first pass, we only consider reads with unique alignments. Counts are recorded in the dict unique.
    In the second pass, we consider the remaining reads. Counts are distributed based on the proportion of uniquely mapped reads.
    '''

    global lower, upper
    lower = frag_len_mean - 3* frag_len_std
    upper = frag_len_mean + 3* frag_len_std

    counts_dict = {} #key = peptide ID, value = object of class Counts
    
    #initiate counts_dict
    for seq_record in SeqIO.parse(reference,"fasta"):
        counts_dict[seq_record.id] = Counts(len(seq_record))

    #first pass
    unique_pass(input_alns,counts_dict)
    
    # print("after first pass")
    # for k,v in counts_dict.items():
    #     print(k)
    #     print(v)

 
    #second pass
    rescue_pass(input_alns,counts_dict)
    # print("after rescue pass")
    # for k,v in counts_dict.items():
    #     print(k)
    #     print(v)

 
    
    #compute TPM
    #read/length
    scaling_factor = 0.0
    for counts in counts_dict.values():
        counts.final_count_norm = counts.final_count/counts.length
        scaling_factor += counts.final_count_norm
    
    # print(scaling_factor)
    
    # print("after normalized final")
    # for k,v in counts_dict.items():
    #     print(k)
    #     print(v)
    
    for counts in counts_dict.values():
        if scaling_factor != 0:
            counts.tpm = counts.final_count_norm * 1000000/scaling_factor

    # print("after tpm")
    # for k,v in counts_dict.items():
    #     print(k)
    #     print(v)
          
    with open(out_counts,"w") as tab_file:
        writer = csv.writer(tab_file, delimiter='\t')

        #adding Header
        writer.writerow(["Name","Length","EffectiveLength","TPM", "NumReads"])

        for key, counts in sorted(counts_dict.items()):
            writer.writerow([key,counts.length,counts.length,counts.tpm,counts.final_count])

    return counts_dict

#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',help="alignments in tab format")
    parser.add_argument('output_file',help="output file containing counts")
    #parser.add_argument('--multi_mapping',action='store_true')
    parser.add_argument('--frag_len_mean',type=float,help="mean of fragment length, when translated. Use the same value as last-pair-probs")
    parser.add_argument('--frag_len_std',type=float,help="standard deviation of fragment length, when translated. Use the same value as last-pair-probs") 
    parser.add_argument('--reference',help="fasta file containing the reference protein sequences")
    
    args = parser.parse_args()
    main(args.input_file,args.output_file, args.frag_len_mean, args.frag_len_std, args.reference)
