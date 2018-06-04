#!/bin/bash

configfile: "config.yaml"

rule all:
     input: "results/Homo_sapiens.GRCh38.92.bed", "results/Homo_sapiens.GRCh38.92.genePred", "data/Homo_sapiens.GRCh38.92.gtf", "results/hg38.chrom.sizes", "exe/gtfToGenePred", "exe/genePredToBed" #expand("results/{accession}.bam", accession=config['all_runs']), 
            # "data/GRCh38","results/multiqc_report.html", expand("results/{accession}_fastqc.html", accession=config['all_runs'])

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
        "results/{accession}.sorted.bam"
    params:
        "-m 4G"
    threads: 8
    wrapper:
        "0.23.1/bio/samtools/sort"

rule download_fetchChromSizes:
     output: "exe/fetchChromSizes"
     shell: "wget --directory-prefix=exe/ http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes"

rule download_gtfToGenePred:
     output: "exe/gtfToGenePred"
     shell: "wget --directory-prefix=exe/ http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred"

rule download_genePredToBed:
     output: "exe/genePredToBed"
     shell: "wget --directory-prefix=exe/ http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed"

rule get_gtf_file:
     output: "data/Homo_sapiens.GRCh38.92.gtf"
     shell: "wget --directory-prefix=data/ ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz && gunzip data/Homo_sapiens.GRCh38.92.gtf.gz"

# derived from https://gist.github.com/gireeshkbogu/f478ad8495dca56545746cd391615b93

rule convert_gtf_to_genepred:
     input:
         script="exe/gtfToGenePred",
         gtf="data/Homo_sapiens.GRCh38.92.gtf"
     output:
         "results/Homo_sapiens.GRCh38.92.genePred"
     shell:
         "{input.script} {input.gtf} {output}"

rule convert_genepred_to_bed12:
     input:
         script="exe/genePredToBed",
         genepred="results/Homo_sapiens.GRCh38.92.genePred"
     output:
         "results/Homo_sapiens.GRCh38.92.bed"
     shell:
         "{input.script} {input.genepred} {output}"

rule get_chrom_sizes:
     input:
         "exe/fetchChromSizes"
     output:
         "results/hg38.chrom.sizes"
     shell:
         "{input} hg38 > results/hg38.chrom.sizes"

rule get_bedgraph:
    input:
        sorted_bam="results/{accession}.sorted.bam",
        genome_sizes="results/hg28.chrom.sizes"
    output:
        "results/{accession}.bedgraph"
    conda:
        "envs/bedtools.yaml"
    benchmark:
        "benchmarks/genomeCoverageBed_{accession}.txt"
    shell:
        "genomeCoverageBed -bg -ibam {input.sorted_bam} -g {input.genome_sizes} -split > {output}"



