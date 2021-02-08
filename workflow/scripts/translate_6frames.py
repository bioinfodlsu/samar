from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def checklen(seq):
    if (len(seq)%3 == 1):
    	return seq + "NN"
    elif(len(seq)%3 == 2):
    	return seq + "N"
    else:
    	return seq


parser = argparse.ArgumentParser()
parser.add_argument('input_file',help="file containing reads")
parser.add_argument('file_type',help="input file type")
parser.add_argument('output_file',help="output file containing translated sequences in all six frames")
args = parser.parse_args()
with open(args.input_file,"r") as fin,open(args.output_file,"w") as fout:

    for seq_record in SeqIO.parse(fin,args.file_type):
        header = seq_record.id
        source = seq_record.seq
        source_rc = source.reverse_complement()

        translated_records = []
        translated_records.append(SeqRecord(seq=checklen(source).translate(), id = header+"_1",description=header+"_1"))
        translated_records.append(SeqRecord(seq=checklen(source[1:]).translate(), id = header+"_2",description=header+"_2"))
        translated_records.append(SeqRecord(seq=checklen(source[2:]).translate(), id = header+"_3",description=header+"_3"))
        translated_records.append(SeqRecord(seq=checklen(source_rc).translate(), id = header+"_4",description=header+"_4"))
        translated_records.append(SeqRecord(seq=checklen(source_rc[1:]).translate(), id = header+"_5",description=header+"_5"))
        translated_records.append(SeqRecord(seq=checklen(source_rc[2:]).translate(), id = header+"_6",description=header+"_6"))
        SeqIO.write(translated_records,fout,"fasta")
