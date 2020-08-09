rule ConvertToFastA:
    input:
        "data/quer/fastq/{sample}.fastq"
    output:
        "data/quer/fasta/{sample}.fasta"
    script:
        "scripts/fastqToFastA.py"
        
rule Align:
    input:
    	quer="data/quer/fasta/{sample}.fasta",
    	ref="data/ref/UP000000803_7227.fasta"
    output:
    	"alignments/{sample}_alg.maf"
    shell:
    	"lastdb flyDB {input.ref} -p \n"
    	"lastal flyDB {input.quer} -p BL62 -F 15 > {output}"
    	
rule ConvertToCSV:
    input:
    	"alignments/{sample}_alg.maf"
    output:
    	"CSV/{sample}.csv"
    script:
    	"scripts/MAFtoCSV.py"
