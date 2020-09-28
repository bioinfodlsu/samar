configfile: "config_file.yaml"
workdir: config["out_dir"]

fasta_interleave = srcdir("scripts/fasta-interleave.sh")
seq_count = srcdir("scripts/seq_count.py")
multimapping: config["multi-mapping"]
def get_read1(wildcards):
    return config["samples"][wildcards.sample_id][0]

def get_read2(wildcards):
    return config["samples"][wildcards.sample_id][1]

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
    shell:
       "lastdb -p {params.index_basename} {input.reference}"

rule align_last:
    input:
        reference_flag = "last_index/index.done",
        reads1 = get_read1,
        reads2 = get_read2
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"])

    output:
        "last_alignments/{sample_id}.tab"

    shell:
        #"{0} {{input.reads1}} {{input.reads2}} | lastal -i1 -p BL62 -F15 {params.last_index_basename} | maf-convert tab > {{output}}".format(os.path.join(workflow.basedir,"scripts2","my-fasta-interleave.sh"))
        "{fasta_interleave}"+" {input.reads1} {input.reads2} | lastal -i1 -p BL62 -F15 {params.last_index_basename} | maf-convert tab > {output}"

rule count_last:
    input:
        alignments = "last_alignments/{sample_id}.tab"
    output:
        "last_counts/{sample_id}/{sample_id}.counts"
    shell:
        "python {seq_count} {input} {output} {multimapping}"
