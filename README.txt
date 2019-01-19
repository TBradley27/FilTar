FilTar - Command Line Tool

Description:

This is a command-line instance of the FilTar tool. Its primary purpose is to provide tissue specific miRNA target prediction in animals. This is achieved by
integrating RNA-Seq data and existing, publicly available algorithms for miRNA target prediction - and allowing users to filter target predictions on the
basis of expression of mRNA transcripts within a given biological context (e.g. tissue or cell line).

It can be executed in two primary modes:

a) Reannotation mode: The RNA-Seq data will be used to reannotate the 3'UTRs for a given biological context (recommended)
b) No reannotation mode: Target predictions will be performed using standard Ensembl/GENCODE annotations

The tool is implemented using the snakemake workflow manager (Koester and Rahmann, 2012 - https://doi.org/10.1093/bioinformatics/bts480) giving it a number of beneficial properities: Extensive automation, modularity, reproducibility, extensibility, configurability, and ease of workflow interpretation and maintenance.

If the tool is executed without implementing an expression threshold or reannotating 3'UTRs then it can function as a highly-automated wrapper for existing
target predictions algorithms such as TargetScan7 (Agarwal et al. 2015 - 10.7554/eLife.05005).

Dependencies:

1. Conda: The tool is a snakemake project and it is therefore essential that the conda package manager is installed on your system. The easiest way to do this is to obtain conda through the relatively lightweight miniconda python distribution, which only contains conda and its dependencies (Installation instructions: https://conda.io/docs/user-guide/install/index.html).

2. Snakemake: Snakemake can be installed using conda essily with the following command (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):

conda install -c bioconda -c conda-forge snakemake

3. filtar (R package): The filtar R package also exists to perform a series of data manipulations within the main FilTar workflow. It can be installed via use of the devtools R package with the following command within the R console:

devtools::install_github('TBradley27/filtar')

4. perl modules: It is not recommended to install perls module using conda*. The recommended method of installing perl modules for us with filtar is to install cpanm perl package manager, and then to manually install perl modules using cpanm. cpanm can itself be installed using conda:

conda install -c bioconda perl-app-cpanminus 

Two perl modules are needed to correctly execute targetscan7 - Bio::TreeIO and Statistics::Lite. They can be installed as follows:

cpanm Bio::Perl
cpanm Statistics::Lite

* This is because conda installs perl modules in different non-interacting directories using different installations of perl, making them very difficult if not possible to use together.

5. Other: For the remaining dependencies, the user has the choice of the following two options:

a)  To allow dependencies to be managed within the snakemake workflow. If the '--use-conda' flag (see later) is used when executing snakemake, then snakemake will use conda to download and/or activate environments containing dependencies needed to execute a given rule (i.e. job within the larger workflow). That environment will then be reactivated whenever that same rule is executed.

For example, is the '--use-conda' flag is used, filtar will download and install HISAT2 and run withins its own self-contained conda environment. The dependency would be installed within conda, so root priviliges would not be needed.

b) The user can manage their own dependencies outside of the filtar workflow. If the '--use-conda' flag is not used, then filtar will try and find HISAT2 within the current environment using the PATH environment variable.

Note 1: Option a) may cause problems on systems such as some HPC environments in which data download is uncoupled from intensive data processing. In which case option b) may be preferable

Note 2: For more flexibility if users would like filtar to manage some dependencies internally (e.g. HISAT2) and some dependencies externally (e.g. trim galore) they could edit the relevant Snakefile and comment out the 'conda' directive from the relevant rule.
In this case, that rule's dependencies would not be managed internally by filtar even when the '--use-conda' flag is used when executing the 'snakemake' command

Usage:

filtar is used as followed:

snakemake {name_of_target_file}

When executing this command, the user must ensure that they are in the project root directory containing the main project makefile.

The name structure of the target file will depend on whether the user wishes to reannotate 3'UTRs or not. For reannotation:

snakemake results/targets/{species}_{tissue}_msa.contextpp.tsv 

The 'species' and 'tissue' wildcards will evaluate to the the species and biological context of interest. The three letter prefix of the species using the first letter of the genus name, and the first two letters of the species name must be used e.g.:

snakemake results/targets/hsa_liver_msa.contextpp.tsv

If the user does not want to reannotate 3'UTR, but still wishes to produce target predictions, then they can run the following command:

snakemake results/targets/no_reannotation/{species}_msa.contextpp.tsv

For targetscan, if the user wishes to identify miRNA target sites, without computing a context++ score (Agarwal et al. 2015 - 10.7554/eLife.05005) for each target site, then the substring 'contextpp' in the target file name should be substituted for 'sites' like so:

snakemake results/targets/no_reannotation/{species}_msa.sites.tsv

If the user wishes to filter computed contextpp scores, then they should add the suffix 'filtered.tsv' to the ususal file name like so:

snakemake results/targets/no_reannotation/{species}_msa.sites.filtered.tsv

Many options can be passed to the snakemake command in order modulate behaviour, reference should be made to the official Snakemake documentation as an exhaustive reference (https://snakemake.readthedocs.io/en/stable/)

Of particular importance, is to note that Snakemake has its own built-in scheduling to manage the execution of different rules. Many rules can be executed in parallel using the '--cores {num_cores}' option.
Combined use of this option and execution within high-performance computing environments enable the execution of rules across many different cores.

Configuration:

The primary configuration which must take place is the following:

a) A mapping between tissue or cell line identifiers, their respective sample accession (e.g. https://www.ncbi.nlm.nih.gov/biosample/), and the mapping of sample accessions to sequencing run accessions
b) A labelling of each run accession to determine if that library corresponds to single-end or paired-end sequencing
c) A selection of which miRNAs (using canonical miRNA names found in miRBase with the three letter species prefix) to use for target prediction. If this configuration is left blank, target prediction is performed for all miRBase annotated miRNAs of that species
d) The expression threshold on which context++ scores are filtered. To manage this, users should edit the 'params' directive of the results/no_reannotation or results/target_predication/targetscan Snakefiles for 'reannotation' and 'no reannotation' analyses respectively

These details can be manually configured using the config/accession_mappings.yaml and config/config.yaml config files respectively. 

For run accessions, users can provide their own data, or aim to download raw data from public repositiories. If the former option is used, then users must manually place their data in either the data/single-end
or the data/paired-end directories depending on which is appropriate. Accessions used do not necessarily need to be of external significance (i.e. universally recognised and understood by other scientists), however if the aim is to download raw sequencing data from public repositories, then run accessions must correspond to the accessions of the run sequence data which is aimed to be downloaded. If user data is used, files names must be in the form {accession}.fastq.gz for single-end data and {accession}_1.fastq.gz and {accession}_2.fastq.gz containing the first and second mate pair reads, respectively.

For the most part, external tools are executed using default parameters which are documented in their respective snakefiles (found within subdirectories of the sub_snakemake directory). Values of parameters can be altered by directly editing the 'param' key-value pairs contained within the rules of these snakefiles.

# WARNING

The script generating the context++ scores uses RNAplfold as a dependency which itself generates ~100KB of data per transcript per species per tissue.
Therefore, generating context++ score for a large number of tissues will consume a lot of disk space. RNAplfold can de deleted after each analysis run
to help mitigate against this issue. ALso, be cautious to delete files if running the same analysis (e.g. same species and tissue) using different 
3'UTR annotations, as targetscan7 will use RNAplfold generated from a previous analysis automatically. Again, manually deleting RNAplfold output will
protect against this issue.














