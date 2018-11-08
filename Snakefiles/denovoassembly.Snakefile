"""----------------------------
Mostafa.Abdel-Glil@fli.de, date: October, 15, 2018, Friedrich-Loeffler-Institut (https://www.fli.de/)
-------------------------------
The pipeline includes the following
- fastqc - quality assessment of raw reads
- sickle - quality trimming
- SPAdes - de novo assembly
- (mini)kraken  - taxonomic classification of fastq reads
                - taxonomic classification of assembled contigs
- Coverage Stats for mapped reads
- quast - quality screen of assemblies
- multiqc report
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
- IMPORTANT: change all run in all rules to shell, so you can execute snakemake with --use-conda
- rule dataset_link, #make sure that the symbolic links are not broken (e.g., if the path is not fully written in config), otherwise copy the files #under test
- mini_kraken and mini_kraken_summary folders, edit the scripts to make changing the folders names felxible
- add log files to the rules
- add benchmark to the rules
- add conda env to the rules
#addrules
- rule theoritical coverage (no * length of reads / assembly size)
-------------------------------
"""
"""
Rules output files
"""
dataset_link=expand(link_dir + "{sample}_1.fastq.gz", sample=sample),
GUNZIP=expand(fastq_dir + "{sample}_1.fastq", sample=sample)
data_preparation =  [dataset_link, GUNZIP]
FastQC=expand(fastqc_dir + "{sample}/{sample}_2_fastqc.zip", sample=sample)
SICKLE=expand(sickle_dir + "{sample}/{sample}_trim_1.fastq.gz", sample=sample)
SPAdes=expand( SPAdes_dir + "{sample}/{sample}.fasta", sample=sample)
Kraken_fastq=expand(kraken_fastq_dir + "{sample}.sequences.kraken", sample=sample)
Kraken_contig=expand( kraken_contig_dir + "{sample}/{sample}.contigs.labels", sample=sample)
CoverageInfo=expand( coverage_info_dir + "{sample}.perBaseCoverage", sample=sample)
Quast=expand( quast_dir + "{sample}/report.html", sample=sample)
denovo_assembly_report= multiqc_dir + "denovo_assembly_report.html"

"""
Raw folder contin a symbolic links for the raw reads
"""
rule dataset_link: #this is to simplify the names of the raw data, so all samples  will be named as sample ID + _1 or _2 for fw and rv reads respect.
    input:
        fw= lambda wildcards: config["samples"][wildcards.sample]["fw"], #fw= lambda wildcards: expand(config["samples"][wildcards.sample]["fw"], sample=sample),
        rv= lambda wildcards: config["samples"][wildcards.sample]["rv"],
    output:
        raw_1=temp(link_dir + "{sample}_1.fastq.gz"),
        raw_2=temp(link_dir + "{sample}_2.fastq.gz")
    run:
        commands = [
            #"""ln -s {input.fw} {output.raw_1} 2> /dev/null && ln -s {input.rv} {output.raw_2} 2> /dev/null || cp -r  {input.fw} {output.raw_1} & cp -r {input.rv} {output.raw_2}"""
            """if  ln -s {input.fw} {output.raw_1} 2> /dev/null && ln -s {input.rv} {output.raw_2} 2> /dev/null ; then \
                echo "sympolic links created"; \
                else echo "sympolic links are not supported. Copying the data...." & cp -r  {input.fw} {output.raw_1} && cp -r {input.rv} {output.raw_2} ; fi"""
            #if the creation of symbolic links is not supported then copy the files
            #To-Do: make sure that symbolic links are not broken, otherwise copy the files #under test
        ]
        for c in commands:
            shell(c)
"""
Decompress the zipped fastq files
"""
rule GUNZIP:
    input:
        fw= lambda wildcards: config["samples"][wildcards.sample]["fw"],
        rv= lambda wildcards: config["samples"][wildcards.sample]["rv"]
    output:
        FASTQ_1= temp(fastq_dir + "{sample}_1.fastq"),
        FASTQ_2= temp(fastq_dir + "{sample}_2.fastq")
    shell:
        "gunzip -c {input.fw} > {output.FASTQ_1} & gunzip -c {input.rv} > {output.FASTQ_2}"

"""
Quality check of raw reads using FastQC
"""
rule FastQC:
    input:
        r1 = link_dir + "{sample}_1.fastq.gz",
        r2 = link_dir + "{sample}_2.fastq.gz"
    output:
        fastqc_zip = fastqc_dir + "{sample}/{sample}_2_fastqc.zip"
    threads: 32
    benchmark:
        benchmarks_folder + "{sample}/fastqc.txt"
    log:
        log_folder + "{sample}/fastqc.log"
    conda:
      envs_folder + "fastqc.yaml"
    run:
        shell("fastqc -t {threads} -f fastq --quiet -o {fastqc_dir}{wildcards.sample} {input.r1} {input.r2} &> {log}")
"""
Quality trimming of raw reads using Sickle
"""
rule Trimming:
    input:
        r1 = link_dir + "{sample}_1.fastq.gz",
        r2 = link_dir + "{sample}_2.fastq.gz"
    output:
        SICKLE_R1= sickle_dir + "{sample}/{sample}_trim_1.fastq.gz",
        SICKLE_R2= sickle_dir + "{sample}/{sample}_trim_2.fastq.gz",
        SICKLE_single= sickle_dir + "{sample}/{sample}_singlesfile.fastq.gz"
    threads: 32
    benchmark:
        benchmarks_folder + "{sample}/sickle.txt"
    log:
        log_folder + "{sample}/sickle.log"
    #conda:
    #  envs_folder + "sickle.yaml"
    shell:
        "sickle pe -l 40 -q 20 -g -t sanger  -f {input.r1} -r {input.r2}  -o {output.SICKLE_R1} -p {output.SICKLE_R2} -s {output.SICKLE_single} &> {log}"
"""
Assembly using SPAdes (`--careful --cov-cutoff auto`)
"""
rule SPAdes: #do assembly with SPAdes,
            #filter contigs for 3x Kmer coverage and 500 bp length, (note: the SPAdes coverage is a kmer coverage, not the real read coverage of the contigs)
            #remove unnecessary folders,
            #collect all filtered contigs from SPAdes into a single folder
    input:
        scripts_dir=config["tools"]["scripts_dir"],
        SICKLE_R1= sickle_dir + "{sample}/{sample}_trim_1.fastq.gz",
        SICKLE_R2= sickle_dir+ "{sample}/{sample}_trim_2.fastq.gz",
        SICKLE_single= sickle_dir + "{sample}/{sample}_singlesfile.fastq.gz"
    output:
        SPAdes_contig=SPAdes_dir + "{sample}/contigs.fasta",
        SPAdes_scaffolds= SPAdes_dir + "{sample}/scaffolds.fasta",
        SPAdes_Filtered_contig= SPAdes_dir + "{sample}/{sample}.fasta",
        collect_contig=FilteredContigs_dir + "{sample}.fasta"
    threads: 32
    benchmark:
        benchmarks_folder + "{sample}/spades.txt"
    log:
        log_folder + "{sample}/spades.log"
    conda:
      envs_folder + "spades.yaml"
    run:
        commands = [
            "spades.py --careful --cov-cutoff auto -t {threads} -1 {input.SICKLE_R1} -2 {input.SICKLE_R2} -s {input.SICKLE_single} -o {SPAdes_dir}{wildcards.sample} &> {log}",
            "python {input.scripts_dir}/filter_contigs.py -i {output.SPAdes_contig} -o  {output.SPAdes_Filtered_contig} -c 3 -l 500",
            "rm -R {SPAdes_dir}{wildcards.sample}/mismatch_corrector && rm -R {SPAdes_dir}{wildcards.sample}/K* && rm -R {SPAdes_dir}{wildcards.sample}/misc",
            "bash {input.scripts_dir}/fa_rename.sh {output.SPAdes_Filtered_contig} {output.collect_contig}"
        ]
        for c in commands:
            shell(c)
"""
Kraken taxonomic assignment of raw reads, screening with mini-kraken DB
"""
rule kraken: #screen the fastq files against minikraken database and summarize the results
    input:
        FASTQ_1=fastq_dir + "{sample}_1.fastq",
        FASTQ_2=fastq_dir + "{sample}_2.fastq",
        db=config["DB"]["minikraken"],
        scripts_dir=config["tools"]["scripts_dir"]
    output:
        mini_kraken_seq= kraken_fastq_dir + "{sample}.sequences.kraken",
        mini_kraken_seq_label= kraken_fastq_dir + "{sample}.sequences.labels",
        mini_kraken_report= kraken_fastq_dir + "{sample}.report.txt",
        mini_kraken_report_summary= kraken_fastq_dir + "{sample}.report_summary.txt"
    threads: 32
    benchmark:
        benchmarks_folder + "{sample}/kraken.txt"
    log:
        log_folder + "{sample}/kraken.log"
    conda:
      envs_folder + "kraken.yaml"
    run:
        commands = [
            "kraken --db {input.db} --paired --check-names --threads {threads} --fastq-input {input.FASTQ_1} {input.FASTQ_2} \
            --output {output.mini_kraken_seq} &> {log}",
            "kraken-translate --db {input.db} {output.mini_kraken_seq} > {output.mini_kraken_seq_label}",
            "kraken-report --db {input.db} {output.mini_kraken_seq} > {output.mini_kraken_report}",
            """ awk '{{print $4}}' {output.mini_kraken_seq_label} | sort | uniq -c | sort -nr > {output.mini_kraken_report_summary}""",
            #see: https://github.com/BacterialCommunitiesAndPopulation/Wednesday18thMay/blob/master/Assembly_Tutorial.md,
            "{input.scripts_dir}/run_kraken_map.sh {kraken_fastq_dir} mini_kraken_summary"
        ]
        for c in commands:
            shell(c)
"""
kraken taxonomic assignment of assembled contigs, screening with kraken DB
see: https://github.com/BacterialCommunitiesAndPopulation/Wednesday18thMay/blob/master/Assembly_Tutorial.md
"""
rule kraken_contig:
    input:
        SPAdes_Filtered_contig= SPAdes_dir + "{sample}/{sample}.fasta",
        db=config["DB"]["kraken"],
        db_krona=config["DB"]["krona"]
    output:
        kraken_contig= kraken_contig_dir + "{sample}/{sample}.contigs.kraken",
        kraken_contig_label= kraken_contig_dir + "{sample}/{sample}.contigs.labels",
        kraken_contig_report_summary= kraken_contig_dir + "{sample}/{sample}.report_summary.txt",
        krona_chart= kraken_contig_dir + "{sample}/{sample}.krona_chart.html"
    threads: 32
    benchmark:
        benchmarks_folder + "{sample}/kraken_contig.txt"
    log:
        log_folder + "{sample}/kraken_contig.log"
    conda:
      envs_folder + "kraken.yaml" #, "krona.yaml"
    run:
        commands = [
            "kraken --db {input.db} --threads {threads} --fasta-input {input.SPAdes_Filtered_contig} --output {output.kraken_contig} &> {log}",
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
        SICKLE_R1= sickle_dir + "{sample}/{sample}_trim_1.fastq.gz",
        SICKLE_R2= sickle_dir + "{sample}/{sample}_trim_2.fastq.gz",
        SPAdes_Filtered_contig= SPAdes_dir + "{sample}/{sample}.fasta"
    output:
        mapped_reads2fasta_SAM=temp(SPAdes_dir + "{sample}/{sample}.fasta.SAM"),
        mapped_reads2fasta_BAM=temp(SPAdes_dir + "{sample}/{sample}.fasta.bam"),
        mapped_reads2fasta_BAM_sorted=temp(SPAdes_dir + "{sample}/{sample}.fasta.bwa.bam"),
        mapping_stats= coverage_info_dir + "{sample}.stats",
        mapping_perBaseCoverage= coverage_info_dir + "{sample}.perBaseCoverage",
        Summary_AverageCoverage= coverage_info_dir + "{sample}.AverageCoverage.txt",
        Summary_contigCoverage= coverage_info_dir + "{sample}.contigCoverage.txt",
        qualimap_report= qualimap_dir + "{sample}/qualimapReport.html"
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
            "/home/software/qualimap/qualimap_v2.2.1/qualimap bamqc -bam {output.mapped_reads2fasta_BAM_sorted} -outdir {qualimap_dir}{wildcards.sample}"
            #if using putty, be sure xming is running, to avoid qualimap error
        ]
        for c in commands:
            shell(c)
"""
quality assessment of contig assembly using QUAST
"""
rule quast:
    input:
        SPAdes_contig= SPAdes_dir + "{sample}/contigs.fasta",
        #SPAdes_scaffolds="SPAdes/{sample}/scaffolds.fasta",
        mapped_reads2fasta_BAM_sorted=temp(SPAdes_dir + "{sample}/{sample}.fasta.bwa.bam")
    output:
        quast_stats= quast_dir + "{sample}/report.html"
    threads: 32
    run:
        commands = [
            "quast.py -t {threads} --bam {input.mapped_reads2fasta_BAM_sorted} -o {quast_dir}{wildcards.sample} -L {input.SPAdes_contig}"
        ]
        for c in commands:
            shell(c)
"""
Summary report using multiqc (fastqc, kraken, quast, qualimap)
"""
rule multiqc:
    input:
        multiqc_bin=config["tools"]["multiqc_bin"],
        FastQC=expand(fastqc_dir + "{sample}/{sample}_2_fastqc.zip", sample=sample),
        mini_kraken_seq=expand( kraken_fastq_dir + "{sample}.sequences.kraken", sample=sample),
        mapping_perBaseCoverage=expand(coverage_info_dir + "{sample}.perBaseCoverage", sample=sample),
        quast_stats=expand(quast_dir + "{sample}/report.html", sample=sample),
    output:
        denovo_assembly_report= multiqc_dir + "denovo_assembly_report.html",
    run:
        commands = [
            "{input.multiqc_bin}/multiqc --filename denovo_assembly_report --force --outdir {multiqc_dir} {fastqc_dir} mini_kraken_summary {qualimap_dir} {quast_dir}",
        ]
        for c in commands:
            shell(c)
