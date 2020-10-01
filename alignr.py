import argparse
import os

#CONFIG
#outdir
#samples
#reference
#multimapping

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument("-cl","--count_last", help="Input is tab alignment file name")
group.add_argument("-al","--align_last", help="Input is tab alignment file name")
group.add_argument("-ld","--last_db", help="Sets last reference sequence",action="store_true")
group.add_argument("-a","--all", help="Runs all rules")

args=parser.parse_args()

#count_last
if(args.count_last):
    os.system("snakemake last_counts/{0}/{0}.counts".format(args.count_last))
#align_last    
elif(args.align_last):
    os.system("snakemake last_alignments/{0}.tab".format(args.align_last)) 
#last_db
elif(args.last_db):
    os.system("snakemake last_index/index.done")
#all
elif(args.all):
    os.system("snakemake")









