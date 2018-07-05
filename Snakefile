#!/bin/bash

configfile: "config.yaml"

rule all:
     input: expand("results/{accession}_trimmed.fq.gz", accession=config['all_brain_runs'])

rule download_sequences:
       input:
       output: 
            "data/{accession, [^_]}.fastq.gz"
       shell: 
            "wget -nv --directory-prefix=data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR{wildcards.accession[3]}{wildcards.accession[4]}{wildcards.accession[5]}/00{wildcards.accession[9]}/{wildcards.accession}/{wildcards.accession}.fastq.gz"

rule download_first_paired_end_file:
       input:
       output:
           "data/{accession}_1.fastq.gz", 
       shell:
          "wget -nv --directory-prefix=data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR{wildcards.accession[3]}{wildcards.accession[4]}{wildcards.accession[5]}/00{wildcards.accession[9]}/{wildcards.accession}/{wildcards.accession}_1.fastq.gz"

rule download_second_paired_end_file:
       input:
       output:
           "data/{accession}_2.fastq.gz",
       shell:
          "wget -nv --directory-prefix=data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR{wildcards.accession[3]}{wildcards.accession[4]}{wildcards.accession[5]}/00{wildcards.accession[9]}/{wildcards.accession}/{wildcards.accession}_2.fastq.gz"

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

rule trim_paired_end_reads:
       input:
          file1="data/{accession}_1.fastq.gz",
          file2="data/{accession}_2.fastq.gz"
       output:
           "results/{accession}_trimmed.fq.gz",
           "results/{accession}.fastq.gz_trimming_report.txt"
       conda:
          "envs/trim-galore.yaml"
       benchmark:
          "benchmarks/trim-galore_{accession}.log"
       shell:
          "trim_galore --output_dir results/ --length 35 --stringency 4 --paired {input.file1} {input.file2}"


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
        threads: 1
        params: 
            index="data/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index",
            extra=""
        wrapper:
            "0.23.1/bio/bowtie2/align"

rule samtools_sort:
    input:
        "results/{accession}.bam"
    output:
        "results/{accession}.bam.sorted"
    params:
        "-m 4G"
    threads: 8
    benchmark:
        "benchmarks/samtools_sort_{accession}.log"
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
     shell: "wget --directory-prefix=data/ ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz && gunzip data/Homo_sapiens.GRCh38.92.gtf.gz && sed -i 's/^/chr/g' data/Homo_sapiens.GRCh38.92.gtf"

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
        sorted_bam="results/{accession}.bam.sorted",
        genome_sizes="results/hg38.chrom.sizes"
    output:
        "results/{accession}.bedgraph"
    conda:
        "envs/bedtools.yaml"
    benchmark:
        "benchmarks/genomeCoverageBed_{accession}.txt"
    shell:
        "genomeCoverageBed -bg -ibam {input.sorted_bam} -g {input.genome_sizes} -split > {output}"

rule extend_3utr_HeLa:
    input:
       script="exe/identifyDistal3UTR.pl",
       bed="results/Homo_sapiens.GRCh38.92.bed",
       bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['HeLa']['mock'])
    output:
       "results/HeLa.utr.bed"
    shell: 
       "{input.script} -i {input.bedgraphs} -m {input.bed} -o {output}"

rule extend_3utrs_HEK293:
    input:
       script="exe/identifyDistal3UTR.pl",
       bed="results/Homo_sapiens.GRCh38.92.bed",
       bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['HEK293']['mock'])
    output:
       "results/HEK293.utr.bed"
    shell:
       "{input.script} -i {input.bedgraphs} -m {input.bed} -o {output}"

rule extend_3utrs_Huh7:
    input:
       script="exe/identifyDistal3UTR.pl",
       bed="results/Homo_sapiens.GRCh38.92.bed",
       bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['Huh7']['mock'])
    output:
       "results/Huh7.utr.bed"
    shell:
       "{input.script} -i {input.bedgraphs} -m {input.bed} -o {output}"

rule extend_3utrs_IMR90:
    input:
       script="exe/identifyDistal3UTR.pl",
       bed="results/Homo_sapiens.GRCh38.92.bed",
       bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['IMR90']['mock'])
    output:
       "results/IMR90.utr.bed"
    shell:
       "{input.script} -i {input.bedgraphs} -m {input.bed} -o {output}"

rule extend_3utrs_brain:
    input:
       script="exe/identifyDistal3UTR.pl",
       bed="results/Homo_sapiens.GRCh38.92.bed",
       bedgraphs= lambda wildcards: 'results/' + str(config['brain'][wildcards.brain_region]) + '.bedgraph'
    output:
       "results/{brain_region}.utr.bed"
    shell:
       "{input.script} -i {input.bedgraphs} -m {input.bed} -o {output}"

rule identify_APA_sites_HeLa:
     input:
        script="exe/predictAPA.pl",
        bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['HeLa']['mock']),
        bed="results/HeLa.utr.bed"
     output:
        "HeLa.APA.txt"
     shell:
        "{input.script} -i {input.bedgraphs} -g 1 -n 2 -u {input.bed}  -o {output}"


rule identify_APA_sites_HEK293:
     input:
        script="exe/predictAPA.pl",
        bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['HEK293']['mock']),
        bed="results/HEK293.utr.bed"
     output:
        "HEK293.APA.txt"
     shell:
        "{input.script} -i {input.bedgraphs} -g 1 -n 2 -u {input.bed}  -o {output}"

rule identify_APA_sites_Huh7:
     input:
        script="exe/predictAPA.pl",
        bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['Huh7']['mock']),
        bed="results/Huh7.utr.bed"
     output:
        "Huh7.APA.txt"
     shell:
        "{input.script} -i {input.bedgraphs} -g 1 -n 4 -u {input.bed}  -o {output}"

rule identify_APA_sites_IMR90:
     input:
        script="exe/predictAPA.pl",
        bedgraphs=expand("results/{accession}.bedgraph", accession=config['PRJNA229375']['IMR90']['mock']),
        bed="results/IMR90.utr.bed"
     output:
        "IMR90.APA.txt"
     shell:
        "{input.script} -i {input.bedgraphs} -g 1 -n 2 -u {input.bed}  -o {output}"

rule identify_APA_sites_brain:
     input:
        script="exe/predictAPA.pl",
        bedgraphs= lambda wildcards: 'results/' + str(config['brain'][wildcards.brain_region]),
        bed="results/{brain_region}.utr.bed"
     output:
        "{brain_region}.APA.txt"
     shell:
        "{input.script} -i {input.bedgraphs} -g 1 -n 1 -u {input.bed}  -o {output}"


rule get_bed6_file: #three_prime_utrs only
        input:
            script="exe/gtf_to_bed.sh",
            gtf="data/Homo_sapiens.GRCh38.92.gtf"
        output:
            "results/Homo_sapiens.GRCh38.92.bed6"
        shell:
            "{input.script} {input.gtf} {output}"

#rule split_bed_files:
#       input: "results/Homo_sapiens.GRCh38.92.bed6"
#       output: "results/Homo_sapiens.GRCh38.92.chr{chrom}.bed6"
#       shell: "awk '$1 == {wildcards.chrom}' {input} {output}"

#rule split_bed_files:
#       input: "results/{cell_line}.utr.bed"
#       output: "results/{cell_line}.utr.chr{chrom}.bed"
#       shell: "awk '$1 == {wildcards.chrom}' {input} {output}"

rule get_extended_bed_file:
         input:
            script="exe/extend_bed.R",
            normal_bed="results/Homo_sapiens.GRCh38.92.bed6",
            extended_bed="results/{cell_line}.utr.bed"
         output:
            "results/{cell_line}.utr.full.bed"
         log:
            "{cell_line}.log"
         shell:
             "Rscript {input.script} {input.normal_bed} {input.extended_bed} {output} 2> {log}"

rule reformat_normal_bed_file:
         input:
           script="exe/reformat_normal_bed.R",
           normal_bed="results/Homo_sapiens.GRCh38.92.bed6"
         output: "results/Homo_sapiens.GRCh38.92.reformatted.bed6"
         shell:
             "Rscript {input.script} {input.normal_bed} {output}"

#WARNING: diff has non-standard exit statuses. A non-zero exit status does not mean programme failure

#rule find_bed_differences:
#	input: 
#           normal_bed="results/Homo_sapiens.GRCh38.92.reformatted.bed6",
#           extended_bed="results/{cell_line}.utr.full.bed"
#	output: "{cell_line}_bed_file_differences.txt"
#	shell: "diff -y --suppress-common-lines {input.normal_bed} {input.extended_bed} > {output}"




