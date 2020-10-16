import argparse
import snakemake

#CONFIG
#outdir
#samples
#reference
#multimapping

#parsers
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument("-cl","--count_last", help="Input is tab alignment file name")
group.add_argument("-al","--align_last", help="Input is tab alignment file name")
group.add_argument("-ld","--last_db", help="Sets last reference sequence",action="store_true")
group.add_argument("-a","--all", help="Runs all rules")

#sub-commands
subparsers = parser.add_subparsers(help='sub-command help')
parser_index = subparsers.add_parser('index', help='constructs the index of the database")
parser_index.add_argument('index', action="store_true")

parser_align = subparsers.add_parser('align', help="aligns reads to the database")
parser_align.add_argument('align', action="store_true")

parser_index = subparsers.add_parser('count', help='counting from the alignments')
parser_index.add_argument('count', action="store_true")

args=parser.parse_args()

#count_last
if(args.count_last or args.count):
    #os.system("snakemake last_counts/{0}/{0}.counts".format(args.count_last))
    snakemake.snakemake("Snakemake",targets=["count_last"])
#align_last    
elif(args.align_last or args.align):
    #os.system("snakemake last_alignments/{0}.tab".format(args.align_last)) 
    snakemake.snakemake("Snakemake",targets=["align_last"])
#last_db
elif(args.last_db or args.index):
    #os.system("snakemake last_index/index.done")
    snakemake.snakemake("Snakemake",targets=["last_db"])
#all
elif(args.all):
    #os.system("snakemake")
    snakemake.snakemake("Snakemake",targets=["all"])









