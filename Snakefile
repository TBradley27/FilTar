#!/bin/bash

configfile: "config.yaml"

rule all:
     input:  expand("results/{accession}.bam", accession=config['all_runs']), 
             "data/GRCh38","results/multiqc_report.html", expand("results/{accession}_fastqc.html", accession=config['all_runs'])

rule download_sequences:
       input:
       output: 
            "data/{accession}.fastq.gz"
       shell: 
            "wget -nv --directory-prefix=data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/00{wildcards.accession[9]}/{wildcards.accession}/{wildcards.accession}.fastq.gz"

rule fastqc:
       input: "data/{accession}.fastq.gz"
       output: 
            "results/{accession}_fastqc.html",
            "results/{accession}_fastqc.zip" 
       conda: "envs/fastqc.yaml"
       shell: "fastqc -o results {input}"

rule trim_reads:
       input: 
           "data/{accession}.fastq.gz"
       output: 
           "results/{accession}_trimmed.fq.gz",
           "results/{accession}.fastq.gz_trimming_report.txt"
       conda: 
          "envs/trim-galore.yaml"
       benchmark:
          "benchmarks/trim-galore_{accession}.log"
       shell:
          "trim_galore --output_dir results/  --length 35 --stringency 4 {input}"

rule multiqc:
       input:
           expand("results/{accession}_trimmed.fq.gz", accession=config['all_runs'])
       output: "results/multiqc_report.html"
       shell: "multiqc -f -o results/ results/"

rule get_bowtie_index:
        output: "data/GRCh38"
        shell:
           "wget --directory-prefix=data/ ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz"

rule map_reads:
        input:
           sample=["results/{accession}_trimmed.fq.gz"]
        output:
           "results/{accession}.bam"
        threads: 12
        params: 
            index="data/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index",
            extra=""
        wrapper:
            "0.23.1/bio/bowtie2/align"

rule samtools_sort:
    input:
        "results/{accession}.bam"
    output:
        "results/{sample}.sorted.bam"
    params:
        "-m 4G"
    threads: 8
    wrapper:
        "0.23.1/bio/samtools/sort"
