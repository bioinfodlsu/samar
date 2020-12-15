from pytools.persistent_dict import PersistentDict
storage = PersistentDict("mystorage")
configfile: "config_file.LAST.sprot.rep5.yaml"
workdir: config["out_dir"]

fastq_to_fasta_interleave = srcdir("scripts/fastq-to-fasta-interleave.sh")
seq_count = srcdir("scripts/seq_count.py")

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
    conda:
    	"env/last.yaml"       
    shell:
       #"lastdb -p {params.index_basename} {input.reference}"
       "lastdb {params.index_basename} {input.reference}"

rule last_frag_stat_est:
    input:
        reference_flag = "last_index/index.done",
        reads1 = get_read1,
        reads2 = get_read2
    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"])
    output :
        frag_len_est="last_alignments/{sample_id}.frag_len_est",
    # shell:
    #     #need if/else to check if user is supplying
    #     """
    #     mkfifo last_alignments/{wildcards.sample_id}_1
    #     mkfifo last_alignments/{wildcards.sample_id}_2
    #     head -n 200000 {input.reads1} > last_alignments/{wildcards.sample_id}_1 &
    #     head -n 200000 {input.reads2} > last_alignments/{wildcards.sample_id}_2 &
    #     {fasta_interleave} last_alignments/{wildcards.sample_id}_1 last_alignments/{wildcards.sample_id}_2 | lastal -i1 -p BL62 -F15 {params.last_index_basename} | last-pair-probs -e > {output.frag_len_est}
    #     rm last_alignments/{wildcards.sample_id}_1 last_alignments/{wildcards.sample_id}_2
    #     """
    conda:
    	"env/last.yaml"
    run:
        #need if/else to check if user is supplying the fragment len info
        shell(
        """
        mkfifo last_alignments/{wildcards.sample_id}_1
        mkfifo last_alignments/{wildcards.sample_id}_2
        gzip -cdf {input.reads1} | head -n 400000 > last_alignments/{wildcards.sample_id}_1 &
        gzip -cdf {input.reads2} | head -n 400000 > last_alignments/{wildcards.sample_id}_2 &
        {fastq_to_fasta_interleave} last_alignments/{wildcards.sample_id}_1 last_alignments/{wildcards.sample_id}_2 | lastal -i1 -p BL62 -F15 {params.last_index_basename} | last-pair-probs -e > {output.frag_len_est}
        rm last_alignments/{wildcards.sample_id}_1 last_alignments/{wildcards.sample_id}_2
        """
        )

        with open(output.frag_len_est) as f:
            for line in f:
                if line.startswith("# estimated mean distance"):
                    mean = line.rstrip().split()[-1]
                    print(mean)
                if line.startswith("# estimated standard deviation of distance"):
                    std = line.rstrip().split()[-1]
                    print(std)

        storage.store("{wildcards.sample_id}_mean",mean)
        storage.store("{wildcards.sample_id}_std",std)


rule align_last:
    input:
        reference_flag = "last_index/index.done",
        reads1 = get_read1,
        reads2 = get_read2,
        frag_len_est="last_alignments/{sample_id}.frag_len_est",

    params:
        last_index_basename="{0}/last_index/index".format(config["out_dir"])

    output:
        "last_alignments/{sample_id}.tab"
    conda:
    	"env/last.yaml"
    run:
        #if fragment length info is not provided
        mean = storage.fetch("{wildcards.sample_id}_mean")
        std = storage.fetch("{wildcards.sample_id}_std")
        shell(
        "{fastq_to_fasta_interleave}"+" {input.reads1} {input.reads2} | lastal -i1 -p BL62 -F15 {params.last_index_basename} |"
        "last-pair-probs -f {mean} -s {std} -m 0.95 |"
        "maf-convert tab > {output}"
        )

rule count_last:
    input:
        alignments = "last_alignments/{sample_id}.tab",
        reference = config["reference_p"]
    output:
        "last_counts/{sample_id}/{sample_id}.counts"
    conda:
    	"env/counting.yaml"
    run:
        mean = storage.fetch("{wildcards.sample_id}_mean")
        std = storage.fetch("{wildcards.sample_id}_std")
        shell("python {seq_count} {input.alignments} {output} --frag_len_mean {mean} --frag_len_std {std} --reference {input.reference}")
        #shell("python {seq_count} {input} {output}")
