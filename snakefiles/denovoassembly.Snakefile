"""----------------------------
Mostafa.Abdel-Glil@fli.de, date: October, 15, 2018
-------------------------------
The pipeline includes the following
- fastqc - quality assessment of raw reads
- sickle - quality trimming
- SPAdes - de novo assembly
- (mini)kraken  - taxonomic classification of fastq reads
                - taxonomic classification of assembled contigs
- Coverage Stats for mapped reads
- quast - quality screen of assemblies
-
-------------------------------
The follwoing code is a modified version, which is based on Snakefiles created by Jörg Linde and Eric Zuchantke at the FLI
also based on the follwoing tutorials
https://slowkow.com/notes/snakemake-tutorial/
https://github.com/leipzig/snakemake-example/blob/master/Snakefile
http://pedagogix-tagc.univ-mrs.fr/courses/ABD/practical/snakemake/snake_intro.html
https://molb7621.github.io/workshop/Classes/snakemake-tutorial.html
https://bitbucket.org/cfce/viper
https://github.com/tanaes/snakemake_assemble/blob/master/bin/snakefiles/assemble
-------------------------------
# To-Do
#edit rules
- add log files to the rules
#addrules
- rule theoritical coverage (no * length of reads/ assembly size)
#Snakefiles
- generate predefined folders
#tools
- define which version of software to be used via specifying the full path
-------------------------------
"""
configfile: "assembly.config.yaml"
sample=config["samples"]
db=config["DB"]["kraken"],
db_krona=config["DB"]["krona"]
scripts_dir=config["tools"]["scripts_dir"]

#directories
fastqc_dir="fastqc/"

#Rules output files
dataset_link=expand("0_fastq.gz/{sample}_1.fastq.gz", sample=sample),
GUNZIP=expand("FASTQ/{sample}_1.fastq", sample=sample)
FastQC=expand("fastqc/{sample}/{sample}_2_fastqc.zip", sample=sample)
SICKLE_R1=expand("SICKLE/{sample}/{sample}_trim_1.fastq.gz", sample=sample)
SPAdes_Filtered_contig=expand("SPAdes/{sample}/{sample}.filtered.contigs.fasta", sample=sample)
mini_kraken_seq=expand("mini_kraken/{sample}.sequences.kraken", sample=sample)
kraken_contig_label=expand("kraken_contig/{sample}/{sample}.contigs.labels", sample=sample)
mapping_perBaseCoverage=expand("Coverage_info/{sample}.perBaseCoverage", sample=sample)
quast_stats=expand("quast/{sample}/report.html", sample=sample)
denovo_assembly_report="denovo_assembly_report/denovo_assembly_report.html"

rule all: #unneeded rules could be deactivated by commenting the line (#)
    input:
        dataset_link,
        GUNZIP,
        FastQC,
        SICKLE_R1,
        SPAdes_Filtered_contig,
        mini_kraken_seq,
        kraken_contig_label,
        mapping_perBaseCoverage,
        quast_stats,
        denovo_assembly_report

"""
Raw folder contin a symbolic links for the raw reads
"""
rule dataset_link: #this is to simplify the names of the raw data, so all samples  will be named as sample ID + _1 or _2 for fw and rv reads respect.
    input:
        fw= lambda wildcards: config["samples"][wildcards.sample]["fw"],
        rv= lambda wildcards: config["samples"][wildcards.sample]["rv"]
    output:
        raw_1=temp("0_fastq.gz/{sample}_1.fastq.gz"),
        raw_2=temp("0_fastq.gz/{sample}_2.fastq.gz")
    shell:
        "ln -s {input.fw} {output.raw_1} &\
        ln -s {input.rv} {output.raw_2}"

"""
Decompress the zipped fastq files
"""
rule GUNZIP:
    input:
        fw= lambda wildcards: config["samples"][wildcards.sample]["fw"],
        rv= lambda wildcards: config["samples"][wildcards.sample]["rv"]
    output:
        FASTQ_1="FASTQ/{sample}_1.fastq",
        FASTQ_2="FASTQ/{sample}_2.fastq"
    shell:
        "gunzip -c {input.fw} > {output.FASTQ_1} &\
        gunzip -c {input.rv} > {output.FASTQ_2}"

"""
Quality check of raw reads using FastQC
"""
rule FastQC:
    input:
        r1 = "0_fastq.gz/{sample}_1.fastq.gz",
        r2 = "0_fastq.gz/{sample}_2.fastq.gz"
    output:
        fastqc_zip = "fastqc/{sample}/{sample}_2_fastqc.zip"
    threads: 32
    run:
        shell("fastqc -t {threads} -f fastq --quiet -o fastqc/{wildcards.sample} {input.r1} {input.r2}")
"""
Quality trimming of raw reads using Sickle
"""
rule Trimming:
    input:
        r1 = "0_fastq.gz/{sample}_1.fastq.gz",
        r2 = "0_fastq.gz/{sample}_2.fastq.gz"
    output:
        SICKLE_R1="SICKLE/{sample}/{sample}_trim_1.fastq.gz",
        SICKLE_R2="SICKLE/{sample}/{sample}_trim_2.fastq.gz",
        SICKLE_single="SICKLE/{sample}/{sample}_singlesfile.fastq.gz"
    threads: 32
    shell:
        "sickle pe -l 40 -q 20 -g -t sanger  -f {input.r1} -r {input.r2}  -o {output.SICKLE_R1} -p {output.SICKLE_R2} -s {output.SICKLE_single}"
"""
Assembly using SPAdes (`--careful --cov-cutoff auto`)
"""
rule SPAdes: #do assembly with SPAdes,
            #filter contigs for 3x Kmer coverage and 500 bp length, (note: the SPAdes coverage is a kmer coverage, not the real read coverage of the contigs)
            #remove unnecessary folders,
            #collect all filtered contigs from SPAdes into a single folder
    input:
        scripts_dir=config["tools"]["scripts_dir"],
        SICKLE_R1="SICKLE/{sample}/{sample}_trim_1.fastq.gz",
        SICKLE_R2="SICKLE/{sample}/{sample}_trim_2.fastq.gz",
        SICKLE_single="SICKLE/{sample}/{sample}_singlesfile.fastq.gz"
    output:
        SPAdes_contig="SPAdes/{sample}/contigs.fasta",
        SPAdes_scaffolds="SPAdes/{sample}/scaffolds.fasta",
        SPAdes_Filtered_contig="SPAdes/{sample}/{sample}.filtered.contigs.fasta",
        collect_contig="FilteredContigs/{sample}.fasta"
    threads: 32
    run:
        commands = [
            "spades.py --careful --cov-cutoff auto -t {threads} -1 {input.SICKLE_R1} -2 {input.SICKLE_R2} -s {input.SICKLE_single} -o SPAdes/{wildcards.sample}",
            "python {input.scripts_dir}/filter_contigs.py -i {output.SPAdes_contig} -o  {output.SPAdes_Filtered_contig} -c 3 -l 500",
            "rm -R SPAdes/{wildcards.sample}/mismatch_corrector && rm -R SPAdes/{wildcards.sample}/K* && rm -R SPAdes/{wildcards.sample}/misc",
            "bash {input.scripts_dir}/fa_rename.sh {output.SPAdes_Filtered_contig} {output.collect_contig}"
        ]
        for c in commands:
            shell(c)
"""
Kraken taxonomic assignment of raw reads, screening with mini-kraken DB
"""
rule kraken: #screen the fastq files against minikraken database and summarize the results
    input:
        FASTQ_1="FASTQ/{sample}_1.fastq",
        FASTQ_2="FASTQ/{sample}_2.fastq",
        db=config["DB"]["minikraken"],
        scripts_dir=config["tools"]["scripts_dir"]
    output:
        mini_kraken_seq="mini_kraken/{sample}.sequences.kraken",
        mini_kraken_seq_label="mini_kraken/{sample}.sequences.labels",
        mini_kraken_report="mini_kraken/{sample}.report.txt",
        mini_kraken_report_summary="mini_kraken/{sample}.report_summary.txt"
    threads: 32
    run:
        commands = [
            "kraken --db {input.db} --paired --check-names --threads {threads} --fastq-input {input.FASTQ_1} {input.FASTQ_2} \
            --output {output.mini_kraken_seq}",
            "kraken-translate --db {input.db} {output.mini_kraken_seq} > {output.mini_kraken_seq_label}",
            "kraken-report --db {input.db} {output.mini_kraken_seq} > {output.mini_kraken_report}",
            """ awk '{{print $4}}' {output.mini_kraken_seq_label} | sort | uniq -c | sort -nr > {output.mini_kraken_report_summary}""",
            #see: https://github.com/BacterialCommunitiesAndPopulation/Wednesday18thMay/blob/master/Assembly_Tutorial.md,
            "{input.scripts_dir}/run_kraken_map.sh mini_kraken/ mini_kraken_summary"
        ]
        for c in commands:
            shell(c)
"""
kraken taxonomic assignment of assembled contigs, screening with kraken DB
see: https://github.com/BacterialCommunitiesAndPopulation/Wednesday18thMay/blob/master/Assembly_Tutorial.md
"""
rule kraken_contig:
    input:
        SPAdes_Filtered_contig="SPAdes/{sample}/{sample}.filtered.contigs.fasta",
        db=config["DB"]["kraken"],
        db_krona=config["DB"]["krona"]
    output:
        kraken_contig="kraken_contig/{sample}/{sample}.contigs.kraken",
        kraken_contig_label="kraken_contig/{sample}/{sample}.contigs.labels",
        kraken_contig_report_summary="kraken_contig/{sample}/{sample}.report_summary.txt",
        krona_chart="kraken_contig/{sample}/{sample}.krona_chart.html"
    threads: 32
    run:
        commands = [
            "kraken --db {input.db} --threads {threads} --fasta-input {input.SPAdes_Filtered_contig} --output {output.kraken_contig}",
            "kraken-translate --db {input.db} {output.kraken_contig} > {output.kraken_contig_label}",
            """ awk '{{print $4}}' {output.kraken_contig_label} | sort | uniq -c | sort -nr > {output.kraken_contig_report_summary}""",
            #see: https://github.com/BacterialCommunitiesAndPopulation/Wednesday18thMay/blob/master/Assembly_Tutorial.md
            "ktImportTaxonomy -q 2 -t 3 -tax {input.db_krona} {output.kraken_contig} -o {output.krona_chart}" #ref: snakefile JLinde
        ]
        for c in commands:
            shell(c)
"""
Coverage, calculate coverage stats after reads mapping to the filtered contigs. Contigs were assembled using SPAdes
"""
rule coverage:#map the fastq files to the assembled contigs and calculate the coverage per base
              #calculate the average coverage per assembly and per contigs
              #run qualimap to be used as input for multiqc
    input:
        scripts_dir=config["tools"]["scripts_dir"],
        SICKLE_R1="SICKLE/{sample}/{sample}_trim_1.fastq.gz",
        SICKLE_R2="SICKLE/{sample}/{sample}_trim_2.fastq.gz",
        SPAdes_Filtered_contig="SPAdes/{sample}/{sample}.filtered.contigs.fasta"
    output:
        mapped_reads2fasta_SAM=temp("SPAdes/{sample}/{sample}.filtered.contigs.fasta.SAM"),
        mapped_reads2fasta_BAM=temp("SPAdes/{sample}/{sample}.filtered.contigs.fasta.bam"),
        mapped_reads2fasta_BAM_sorted=temp("SPAdes/{sample}/{sample}.filtered.contigs.fasta.bwa.bam"),
        mapping_stats="Coverage_info/{sample}.stats",
        mapping_perBaseCoverage="Coverage_info/{sample}.perBaseCoverage",
        Summary_AverageCoverage="Coverage_info/{sample}.AverageCoverage.txt",
        Summary_contigCoverage="Coverage_info/{sample}.contigCoverage.txt",
        qualimap_report="qualimap/{sample}/qualimapReport.html"
    threads: 32
    run:
        commands = [
            "bwa index {input.SPAdes_Filtered_contig} ",
            "bwa mem -t {threads} {input.SPAdes_Filtered_contig} {input.SICKLE_R1} {input.SICKLE_R2} > {output.mapped_reads2fasta_SAM}",
            "samtools view -@ {threads} -T {input.SPAdes_Filtered_contig} -bS -o {output.mapped_reads2fasta_BAM} {output.mapped_reads2fasta_SAM}",
            "samtools sort -@ {threads} -T {input.SPAdes_Filtered_contig}.bwa -o {output.mapped_reads2fasta_BAM_sorted} {output.mapped_reads2fasta_BAM}",
            "samtools stats {output.mapped_reads2fasta_BAM_sorted} > {output.mapping_stats}",
            "samtools depth {output.mapped_reads2fasta_BAM_sorted} > {output.mapping_perBaseCoverage}",
            """awk '{{cnt+=$3; n+=1}} END{{ if (n > 0){{ print cnt/n }} else {{ print "0" }} }}' {output.mapping_perBaseCoverage} > {output.Summary_AverageCoverage}""",
            "{input.scripts_dir}/contigCoverage.sh  {output.mapping_perBaseCoverage} > {output.Summary_contigCoverage}",
            "/home/software/qualimap/qualimap_v2.2.1/qualimap bamqc -bam {output.mapped_reads2fasta_BAM_sorted} -outdir qualimap/{wildcards.sample}"
            #if using putty, be sure xming is running, to avoid qualimap error
        ]
        for c in commands:
            shell(c)
"""
quality assessment of contig assembly using QUAST
"""
rule quast:
    input:
        SPAdes_contig="SPAdes/{sample}/contigs.fasta",
        #SPAdes_scaffolds="SPAdes/{sample}/scaffolds.fasta",
        mapped_reads2fasta_BAM_sorted=temp("SPAdes/{sample}/{sample}.filtered.contigs.fasta.bwa.bam")
    output:
        quast_stats="quast/{sample}/report.html"
    threads: 32
    run:
        commands = [
            "quast.py -t {threads} --bam {input.mapped_reads2fasta_BAM_sorted} -o quast/{wildcards.sample} -L {input.SPAdes_contig}"
        ]
        for c in commands:
            shell(c)
"""
Summary report using multiqc (fastqc, kraken, quast, qualimap)
"""
rule multiqc:
    input:
        multiqc_bin=config["tools"]["multiqc_bin"],
        FastQC=expand("fastqc/{sample}/{sample}_2_fastqc.zip", sample=sample),
        mini_kraken_seq=expand("mini_kraken/{sample}.sequences.kraken", sample=sample),
        mapping_perBaseCoverage=expand("Coverage_info/{sample}.perBaseCoverage", sample=sample),
        quast_stats=expand("quast/{sample}/report.html", sample=sample),
    output:
        denovo_assembly_report="denovo_assembly_report/denovo_assembly_report.html",
    run:
        commands = [
            "{input.multiqc_bin}/multiqc --filename denovo_assembly_report --force --outdir denovo_assembly_report fastqc mini_kraken_summary qualimap quast",
        ]
        for c in commands:
            shell(c)
