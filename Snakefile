import glob
import os
import re

configfile: "config.yaml"
CDHIT_THR = config["cdhit"]
CDHIT_THR = 0.8

HITS=[os.path.splitext(os.path.basename(filepath))[0] for filepath in glob.iglob('data/*.hits')]

#print(HITS)

wildcard_constraints:
    CDHIT_THR="[0-9]\.?[0-9]?",
    i="\d+"

rule all:
    input:
        #expand("extracted/{hit}.fasta", hit=HITS)
        #expand("aggregated/{CDHIT_THR}/{hit}.hmm", hit=HITS, CDHIT_THR=CDHIT_THR),
        expand("final_hmms/combined_{CDHIT_THR}.hmm", CDHIT_THR=CDHIT_THR),
        #"aggregated/combined_{{CDHIT_THR}}.hmm"

rule extract:
    input:
        fa=config["extract"]["input_fasta"],
        q="data/{hit}.hits"
    output:
        "extracted/{hit}.fasta"
    shell:
       "seqkit grep --pattern-file {input.q} {input.fa} > {output}"

rule cdhit:
    input:
        "extracted/{hit}.fasta"
    params:
        cores=config["general"]["cores"],
        memory=config["general"]["memory"]
    output:
         cl="cdhit/{CDHIT_THR}/{hit}",
         clstr="cdhit/{CDHIT_THR}/{hit}.clstr"
    shell:
         "cd-hit -i {input} -o {output.cl} -c {CDHIT_THR} -n 5 -T {params.cores} -M {params.memory} -d 0"

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint clustering:
    input:
        clstr="cdhit/{CDHIT_THR}/{hit}.clstr",
        db="extracted/{hit}.fasta"
    output:
        dir=directory("cluster_fasta/{CDHIT_THR}/{hit}"),
        #clusters=directory("clustering/{hit}")
    shell:
        "make_multi_seq.pl {input.db} {input.clstr} {output.dir} 5 #; find {output.dir} -empty -type d -delete"

# an intermediate rule
rule msa:
    input:
        "cluster_fasta/{CDHIT_THR}/{hit}/{i}"
    output:
        "msa_bundled/{CDHIT_THR}-{hit}-{i}.clw"
        #"post/{hit}/{i}.txt"i
    shell:
        #"cp {input} {output}"
        "muscle -in {input} -out {output} -clw"

rule hmm_build:
    input:
        "msa_bundled/{CDHIT_THR}-{hit}-{i}.clw"
    output:
        "hmms/{CDHIT_THR}-{hit}-{i}.hmm"
    shell:
        "hmmbuild {output} {input}"

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
    a = expand("hmms/{CDHIT_THR}-{hit}-{i}.hmm",
    hit=wildcards.hit,
    CDHIT_THR=wildcards.CDHIT_THR,
    i=glob_wildcards(os.path.join(checkpoint_output, "{i}")).i)
    regex = re.compile(r'.*?\d+\.hmm')
    my_list = list(filter(regex.match, a))
    print(my_list)
    if my_list:
        return(my_list)
    else:
        return["./default_file"]

# an aggregation over all produced clusters
rule aggregate:
    input:
        aggregate_input
    output:
       #combined="aggregated/combined_{{CDHIT_THR}}.hmm",
       "aggregated/{CDHIT_THR}/{hit}.hmm"
    shell:
        "cat {input} > {output}"
        #"cat {input} > {output.combined}"

# an aggregation over all produced clusters
rule combine:
    input:
         ["aggregated/{{CDHIT_THR}}/{hit}.hmm".format(hit=hit) for hit in HITS]
    output:
       "final_hmms/combined_{CDHIT_THR}.hmm"
    shell:
        "cat {input} > {output}"
        #"cat {input} > {output.combined}"


rule clean:
    shell:
        '''
        rm -rf aggregated
        '''