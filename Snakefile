configfile: "config_file.LAST.sprot.rep5.yaml"
workdir: config["out_dir"]

fastq_to_fasta_interleave = srcdir("scripts/fastq-to-fasta-interleave.sh")
seq_count = srcdir("scripts/seq_count.py")
trans_6frame = srcdir("scripts/translate_6frames.py")

def get_read1(wildcards):
    return config["samples"][wildcards.sample_id][0]

def get_read2(wildcards):
    return config["samples"][wildcards.sample_id][1]
    
def getMean(wildcard):
    try:
	    with open("last_alignments/{}.frag_len_est".format(wildcard)) as f:
	    	for line in f:
	    		if line.startswith("# estimated mean distance"):
	    			return(line.rstrip().split()[-1])
    except FileNotFoundError:
    	return(80)

def getSTD(wildcard):
    try:
    	with open("last_alignments/{}.frag_len_est".format(wildcard)) as f:
    		for line in f:
    			if line.startswith("# estimated standard deviation of distance"):
		    		return(line.rstrip().split()[-1])
    except FileNotFoundError:
    	return(8)

rule all:
    input:
        [
            expand("last_counts/{sample_id}/{sample_id}.counts",sample_id=config["samples"].keys())
        ]
	
rule last_db:
    input:
        reference = config["reference_p"]
    output:
       touch('last_index/index.done')
    params:
       index_basename="{0}/last_index/index".format(config["out_dir"])
    conda:
    	"env/last.yaml"       
    shell:
       "lastdb -p {params.index_basename} {input.reference}"
       #"lastdb {params.index_basename} {input.reference}"

rule last_scoring:
    input:
    	reference_flag = "last_index/index.done",
	query = config["samples"][list(config["samples"].keys())[0]][0],
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"])
    output:
    	"last_sample/scoring_scheme"
    conda:
    	"env/lastxpy.yaml"
    shell:
    	"""
    	head -n 400000 {input.query} > last_sample/sample_1.sample.fa
    	python {trans_6frame} last_sample/sample_1.sample.fa  last_sample/translated.sample_1.sample.fa
    	last-train --revsym --matsym --gapsym --sample-number=600000 -S0 {params.last_index_basename} last_sample/translated.sample_1.sample.fa> {output}
	"""
rule last_frag_stat_est:
    input:
        reference_flag = "last_index/index.done",
        reads1 = get_read1,
        reads2 = get_read2,
        scoring = "last_sample/scoring_scheme",
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"])
    output :
        frag_len_est="last_alignments/{sample_id}.frag_len_est",
    conda:
    	"env/last.yaml"
    shell:
        """
        mkfifo last_alignments/{wildcards.sample_id}_1
        mkfifo last_alignments/{wildcards.sample_id}_2
        gzip -cdf {input.reads1} | head -n 400000 > last_alignments/{wildcards.sample_id}_1 &
        gzip -cdf {input.reads2} | head -n 400000 > last_alignments/{wildcards.sample_id}_2 &
        {fastq_to_fasta_interleave} last_alignments/{wildcards.sample_id}_1 last_alignments/{wildcards.sample_id}_2 | lastal -i1 -p {input.scoring} -F15 {params.last_index_basename} | last-pair-probs -e > {output.frag_len_est}
        rm last_alignments/{wildcards.sample_id}_1 last_alignments/{wildcards.sample_id}_2
        """

rule align_last:
    input:
        reference_flag = "last_index/index.done",
        reads1 = get_read1,
        reads2 = get_read2,
        frag_len_est="last_alignments/{sample_id}.frag_len_est",
        scoring = "last_sample/scoring_scheme",
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"]),
        mean = getMean,
        std = getSTD,

    output:
        "last_alignments/{sample_id}.tab"
    conda:
    	"env/last.yaml"
    shell:
        "{fastq_to_fasta_interleave}"+" {input.reads1} {input.reads2} | lastal -i1 -p {input.scoring} -F15 {params.last_index_basename} |"
        "last-pair-probs -f {params.mean} -s {params.std} -m 0.95 |"
        "maf-convert tab > {output}"

rule count_last:
    input:
        alignments = "last_alignments/{sample_id}.tab",
        reference = config["reference_p"]
    params:
    	mean = getMean,
    	std = getSTD,
    output:
        "last_counts/{sample_id}/{sample_id}.counts"
    conda:
    	"env/counting.yaml"
    shell:
        "python {seq_count} {input.alignments} {output} --frag_len_mean {params.mean} --frag_len_std {params.std} --reference {input.reference}"
        
