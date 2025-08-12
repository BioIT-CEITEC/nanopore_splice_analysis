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
ANNOTATION = config["annotation"]

rule all:
    input:
        expand("splice_analysis/{sample_name}/{sample_name}_talon_read_annot.tsv", sample_name = sample_tab.sample_name),
        expand("splice_analysis/{sample_name}/{sample_name}_filtered_read_annot.tsv", sample_name = sample_tab.sample_name),
        expand("splice_analysis/{sample_name}/{sample_name}_talon.gtf", sample_name = sample_tab.sample_name),
        expand("splice_analysis/{sample_name}/{sample_name}_talon_abundance.tsv", sample_name = sample_tab.sample_name)


# rule convert_bam_to_sam: 
#     input:expand('aligned/{sample_name}/{sample_name}_sorted.bam', sample_name = sample_tab.sample_name)
#     output: 'aligned/{sample_name}/{sample_name}_sorted.sam'
#     conda: "envs/talon.yaml"
#     shell:
#         """
#         samtools view -h -o {output} {input}
#         """

rule label_reads:
    input: expand("aligned/{sample_name}/{sample_name}_sorted.bam", sample_name = sample_tab.sample_name)
    output: 'labeled/{sample_name}/{sample_name}_labeled.sam'
    params:
        genome = config["organism_fasta"],
        output_dir = lambda wildcards: f"labeled/{wildcards.sample_name}/{wildcards.sample_name}"
    conda: "envs/talon.yaml"
    threads: workflow.cores * 0.75
    shell:
        """
        mkdir -p {params.output_dir}
        talon_label_reads --f={input} --t {threads} --g={params.genome} --deleteTmp  --o={params.output_dir}
        """    

rule initialize_talon_database:
    input:
        ref_gtf = config["organism_gtf"]
    output: "splice_analysis/{sample_name}/talon.db"
    conda: "envs/talon.yaml"
    params:
        annotation=ANNOTATION, 
        genome="hg38",
        output_prefix="splice_analysis/{sample_name}/talon"
    shell:
        """
        talon_initialize_database \
            --f {input.ref_gtf} \
            --g {params.genome} \
            --a  {params.annotation} \
            --o {params.output_prefix} 
        """

rule create_talon_config:
    input:
        sam_files=expand("labeled/{sample_name}/{sample_name}_labeled.sam", sample_name = sample_tab.sample_name)
    output:
        "splice_analysis/config.csv"
    run:
        import os

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as f:
            for file in input.sam_files:
                base = os.path.basename(file).removesuffix(".sam")
                sample = os.path.basename(file).replace("_labeled.sam", "")
                f.write(f"{base},{sample},ONT,{file}\n")

rule talon_annotate:
    input: 
        config="splice_analysis/config.csv",
        db="splice_analysis/{sample_name}/talon.db"
    output: 
        tsv="splice_analysis/{sample_name}/{sample_name}_talon_read_annot.tsv"
    params:
        genome="hg38",
        output_prefix="splice_analysis/{sample_name}/{sample_name}"
    threads: workflow.cores * 0.75
    conda: "envs/talon.yaml"
    shell:
        """
        talon \
            --f {input.config} \
            --db {input.db} \
            --threads {threads}  \
            --build {params.genome} \
            --o {params.output_prefix} \
        """

rule talon_filter_transcripts:
    input:
        db="splice_analysis/{sample_name}/talon.db",
        tsv="splice_analysis/{sample_name}/{sample_name}_talon_read_annot.tsv"
    output:
        filtered="splice_analysis/{sample_name}/{sample_name}_filtered_read_annot.tsv"
    params:
        annotation=ANNOTATION
    conda: "envs/talon.yaml"
    shell:
        """
        talon_filter_transcripts \
            --db {input.db} \
            -a {params.annotation} \
            --o {output.filtered} 
        """

rule talon_create_gtf:
    input:
        db="splice_analysis/{sample_name}/talon.db",
        tsv="splice_analysis/{sample_name}/{sample_name}_talon_read_annot.tsv"
    output:
        "splice_analysis/{sample_name}/{sample_name}_talon.gtf"
    params:
        annotation=ANNOTATION, 
        genome="hg38",
        output_prefix="splice_analysis/{sample_name}/{sample_name}"
    conda: "envs/talon.yaml"
    shell:
        """
        talon_create_GTF --db {input.db} -a {params.annotation} -b {params.genome} --o {params.output_prefix}
        """

rule talon_abundance:
    input:
        db="splice_analysis/{sample_name}/talon.db",
        tsv="splice_analysis/{sample_name}/{sample_name}_talon_read_annot.tsv"
    output:
        "splice_analysis/{sample_name}/{sample_name}_talon_abundance.tsv"
    params:
        annotation=ANNOTATION, 
        genome="hg38",
        output_prefix="splice_analysis/{sample_name}/{sample_name}"
    conda: "envs/talon.yaml"
    shell:
        """
        talon_abundance \
            --db {input.db} \
            -a {params.annotation} \
            -b {params.genome} \
            --o {params.output_prefix} \
        """
