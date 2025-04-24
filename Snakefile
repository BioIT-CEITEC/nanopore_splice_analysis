from pathlib import Path
import pandas as pd

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"] 

##### BioRoot utilities - reference #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*
config = BR.load_organism()
sample_tab = BR.load_sample()

rule all:
    input:
        expand("splice_analysis/{sample_name}/{sample_name}_filtered_read_annot.tsv", sample_name = sample_tab.sample_name),
        expand("splice_analysis/{sample_name}/{sample_name}_transcripts.gtf", sample_name = sample_tab.sample_name),
        expand("splice_analysis/{sample_name}/{sample_name}_abundance.tsv", sample_name = sample_tab.sample_name)


rule convert_bam_to_sam: 
    input:'aligned/{sample_name}/{sample_name}_sorted.bam'
    output: 'aligned/{sample_name}/{sample_name}_sorted.sam'
    conda: "envs/talon.yaml"
    shell:
        """
        samtools view -h -o {output} {input}
        """

rule label_reads:
    input:'aligned/{sample_name}/{sample_name}_sorted.sam'
    output: 'aligned/{sample_name}/labeled_{sample_name}_sorted.sam'
    params:
        genome = config["organism_fasta"]
    conda: "envs/talon.yaml"
    threads: workflow.cores * 0.75
    shell:
        """
        talon_label_reads --f={input} --t {threads} --g={params.genome} --deleteTmp  --o="labeled_"
        """    

rule initialize_talon_database:
    input:
        ref_gtf = config["organism_gtf"]
    output: "splice_analysis/{sample_name}/talon.db"
    conda: "envs/talon.yaml"
    threads: workflow.cores * 0.75
    params:
        annotation="r110", 
        genome="hg38"
    shell:
        """
        talon_initialize_database \
            --f {input.ref_gtf} \
            --g {params.genome} \
            --a  {params.annotation} \
            --l ont \
            --o {output}
        """

rule talon_annotate:
    input: 
        sam='aligned/{sample_name}/labeled_{sample_name}_sorted.sam',
        db="splice_analysis/{sample_name}/talon.db"
    output: 
        tsv="splice_analysis/{sample_name}/{sample_name}read_annot.tsv"
    params:
        sample_name=lambda wildcards: wildcards.sample_name,
        genome="hg38"
    threads: workflow.cores * 0.75
    conda: "envs/talon.yaml"
    shell:
        """
        talon \
            --f {input.sam} \
            --db {input.db} \
            --threads {threads}  \
            --build {params.genome} \
            --o {params.sample_name} \
        """

rule talon_filter_transcripts:
    input:
        db="splice_analysis/{sample_name}/talon.db"
    output:
        filtered="splice_analysis/{sample_name}/{sample_name}_filtered_read_annot.tsv"
    params:
        annotation="r110"
    conda: "envs/talon.yaml"
    shell:
        """
        talon_filter_transcripts \
            --db {input.db} \
            --a {params.annotation} \
            --o {output.filtered} 
        """

rule talon_create_gtf:
    input:
        db="splice_analysis/{sample_name}/talon.db"
    output:
        "splice_analysis/{sample_name}/{sample_name}_transcripts.gtf"
    params:
        annotation="r110"
    conda: "envs/talon.yaml"
    shell:
        """
        talon_create_GTF \
            --db {input.db} \
            --a {params.annotation} \
            --o {output}
        """

rule talon_abundance:
    input:
        db="splice_analysis/{sample_name}/talon.db"
    output:
        "splice_analysis/{sample_name}/{sample_name}_abundance.tsv"
    params:
        annotation="r110",
        genome="hg38"
    conda: "envs/talon.yaml"
    shell:
        """
        talon_abundance \
            --db {input.db} \
            --a {params.annotation} \
            --b {params.genome} \
            --o splice_analysis/{wildcards.sample_name}/{wildcards.sample_name}_abundance
        """
