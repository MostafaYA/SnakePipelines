"""----------------------------
Mostafa.Abdel-Glil@fli.de, date: October, 15, 2018
-------------------------------
The pipeline includes the following
- ariba AMR screnning  - Antimicrobial resistance
- ariba VF screnning - Virulence factors screening
-------------------------------
# To-Do
#addrules
- AMR
- VF
-------------------------------
"""
configfile: "assembly.config.yaml"
sample=config["samples"]

#Rules output files
#aribaAMR=expand( prokka_dir + "{sample}", sample=sample)
#Roary= Roary_dir + "pan_genome_reference.fa"

"""
Annotation using Prokka
"""
rule Prokka: #run prokka for the annotation
            #aggregate the gff files from the prokka output
    input:
        #prokka_bin=config["tools"]["prokka_bin"],
        collect_contig= FilteredContigs_dir + "{sample}.fasta",
    output:
        #prokka_gff= prokka_dir + "{sample}/{sample}.gff",
        prokka_outdir=directory(prokka_dir + "{sample}"),
    threads: 16
    benchmark: "benchmarks/{sample}/prokka.txt"
    log: "Log/{sample}/prokka.log"
    conda: "envs/prokka.yaml" #needs revision
    params:
        #prokka_outdir=directory(prokka_dir + "{sample}"),
        prefix="{sample}",
        strain_name="{sample}",
        locustag= lambda wildcards: config["samples"][wildcards.sample]["locustag"],
        genustag= lambda wildcards: config["samples"][wildcards.sample]["genustag"],
    run:
        commands = [
        "prokka --cpus {threads} --kingdom Bacteria --outdir {output.prokka_outdir}/ --prefix {params.prefix} --genus {params.genustag} \
                --locustag {params.locustag} --strain {params.strain_name} {input.collect_contig} &> {log}",
        "mkdir -p {gff_dir} && ln -sfn $(realpath {output.prokka_outdir}/{wildcards.sample}.gff 2>/dev/null) {gff_dir} ||\
                cp -r {output.prokka_outdir}/{wildcards.sample}.gff {gff_dir}",
        #mk a dir for gff files if the dir does not exist, then create a sympolic link for the gff files. if necessary, force the creation via updating the existed link.
        #If linking doesnt work then copy the data
        ]
        for c in commands:
            shell(c)
