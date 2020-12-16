import argparse
import snakemake
import yaml

#CONFIG
#outdir
#samples
#reference
#multimapping

def yamlInp(inp, val):
    with open("config.yaml") as f:
    	list_doc = yaml.load(f,Loader=yaml.BaseLoader)
    
    for key in list_doc.keys():
    	if(key == inp):
    	    list_doc[inp] = val
    	    
    with open("config.yaml","w") as f:
    	yaml.dump(list_doc,f)

#parsers
parser = argparse.ArgumentParser()
parser.add_argument("-ref","--reference", help="Uses input as reference for pipeline")
parser.add_argument("-red","--reads", help="Uses input as query for pipeline")
parser.add_argument("-out","--out_dir", help="Directory where output will be located")

#parser.add_argument("-ref","--reference", help="Uses input as reference for pipeline")
	
group = parser.add_mutually_exclusive_group()
group.add_argument("-cl","--count_last", help="Counts the amount of reads aligned to each transcript",action="store_true")
group.add_argument("-al","--align_last", help="Aligns the input query and reference sequences",action="store_true")
group.add_argument("-ld","--last_db", help="Sets the reference sequence",action="store_true")
group.add_argument("-a","--all", help="Sets the reference sequence, aligns queries, and counts alignments",action="store_true")
	
'''
#sub-commands
subparsers = parser.add_subparsers(help='sub-command help')
parser_index = subparsers.add_parser('index', help='constructs the index of the database')
parser_index.add_argument('index', action="store_true")

parser_align = subparsers.add_parser('align', help="aligns reads to the database")
parser_align.add_argument('align', action="store_true")

parser_index = subparsers.add_parser('count', help='counting from the alignments')
parser_index.add_argument('count', action="store_true")
'''
args=parser.parse_args()

if(args.reference):
    print("reference")
    yamlInp("reference",args.reference)

if(args.reads):
    print("reads")
    yamlInp("reads",args.reads)

if(args.out_dir):
    print("out")
    yamlInp("out",args.out_dir)

#count_last
if(args.count_last):
    #os.system("snakemake last_counts/{0}/{0}.counts".format(args.count_last))
    snakemake.snakemake("Snakefile",targets=["count_last"],use_conda=True)
    #print("count")
#align_last    
elif(args.align_last):
    #os.system("snakemake last_alignments/{0}.tab".format(args.align_last)) 
    snakemake.snakemake("Snakefile",targets=["align_last"],use_conda=True)
    #print("align")
#last_db
elif(args.last_db):
    #os.system("snakemake last_index/index.done")
    snakemake.snakemake("Snakefile",targets=["last_db"],use_conda=True)
    #print("db")
#all
elif(args.all):
    #os.system("snakemake")
    snakemake.snakemake("Snakefile",targets=["all"],use_conda=True)
    #print("all")

		
    	






