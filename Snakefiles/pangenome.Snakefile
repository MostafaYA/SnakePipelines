"""----------------------------
Mostafa.Abdel-Glil@fli.de, date: October, 15, 2018
-------------------------------
The pipeline includes the following
- Prokka - Annotation
- Roary - Pangenome
-------------------------------
# To-Do
- Prokka, write the path of envs in the config. Accordingly, capture the envs into the rules
#addrules
- Roary
- multiqc report for the annotation
- parse Roary output
-------------------------------
"""
#Rules output files
Prokka=expand( prokka_dir + "{sample}/{sample}.gff", sample=sample)
Roary= Roary_dir + "pan_genome_reference.fa"
#collect_GFF=expand( gff_dir + "{sample}.gff", sample=sample)
#Prokka_addhocAssemblies=expand( prokka_assembly_dir + "{genome}", genome=genome)
"""
Annotation using Prokka
"""
rule Prokka: #run prokka for the annotation
            #aggregate the gff files from the prokka output
    input:
        #prokka_bin=config["tools"]["prokka_bin"],
        collect_contig= FilteredContigs_dir + "{sample}.fasta",
    output:
        prokka_gff= prokka_dir + "{sample}/{sample}.gff",
        #prokka_outdir=directory(prokka_dir + "{sample}"),
    threads: 32
    benchmark: "benchmarks/{sample}/prokka.txt"
    log: "Log/{sample}/prokka.log"
    conda: "envs/prokka.yaml" #needs revision
    params:
        prokka_outdir= prokka_dir + "{sample}",
        prefix="{sample}",
        strain_name="{sample}",
        locustag= lambda wildcards: config["samples"][wildcards.sample]["locustag"],
        genustag= lambda wildcards: config["samples"][wildcards.sample]["genustag"],
    run:
        commands = [
        "prokka --cpus {threads} --force --kingdom Bacteria --outdir {params.prokka_outdir}/ --prefix {params.prefix} --genus {params.genustag} \
                --locustag {params.locustag} --strain {params.strain_name} {input.collect_contig} &> {log}",
        ]
        for c in commands:
            shell(c)
"""
Annotation using Prokka for the already assembled genomes

rule Prokka_addhocAssemblies: #run prokka for the annotation
            #aggregate the gff files from the prokka output
    input:
        #prokka_bin=config["tools"]["prokka_bin"],
        adhoc_assemblies=lambda wildcards: config["genomes"][wildcards.genome]["fasta"]
    output:
        #prokka_gff= prokka_dir + "{sample}/{sample}.gff",
        prokka_outdir_adhoc_assemblies=directory(prokka_dir + "{genome}"),
    threads: 32
    benchmark: "benchmarks/{genome}/prokka.txt"
    log: "Log/{genome}/prokka.log"
    conda: "envs/prokka.yaml" #needs revision
    params:
        prefix="{genome}",
        strain_name="{genome}",
        locustag= lambda wildcards: config["samples"][wildcards.genome]["locustag"],
        genustag= lambda wildcards: config["samples"][wildcards.genome]["genustag"],
    run:
        commands = [
        "prokka --cpus {threads} --kingdom Bacteria --outdir {output.prokka_outdir_adhoc_assemblies}/ --prefix {params.prefix} --genus {params.genustag} \
                --locustag {params.locustag} --strain {params.strain_name} {input.adhoc_assemblies} &> {log}",
        "mkdir -p {gff_dir} && ln -sfn $(realpath {output.prokka_outdir_adhoc_assemblies}/{wildcards.genome}.gff 2>/dev/null) {gff_dir} ||\
                cp -r {output.prokka_outdir_adhoc_assemblies}/{wildcards.genome}.gff {gff_dir}",
        #mk a dir for gff files if the dir does not exist, then create a sympolic link for the gff files. if necessary, force the creation via updating the existed link.
        #If linking doesnt work then copy the data
        ]
        for c in commands:
            shell(c)
"""
"""
Pangenome using Roary
"""
rule Roary: #run roary
    input:
        scripts_dir=config["tools"]["scripts_dir"],
        prokka_gff= expand(prokka_dir + "{sample}/{sample}.gff", sample=sample)
        #gff_dir=directory(gff_dir)
        #expand( gff_dir + "{sample}.gff", sample=sample)
        #prokka_outdir=directory(prokka_dir + "{sample}"),
    output:
        #Roary_outdir=directory(Roary_dir),
        tmp_dir=temp(directory(temporary_todelete)),
        Roary_pangenome_fa= Roary_dir + "pan_genome_reference.fa"
    threads: 64
    benchmark: "benchmarks/roary.txt"
    log: "Log/roary.log"
    conda: "envs/roary.yaml" #needs revision
    #params:
        #Roary_outdir=directory(Roary_dir),
    run:
        commands = [
        "bash {input.scripts_dir}/fixRoaryOutDirError.sh {Roary_dir} {output.tmp_dir}",
        "roary -p {threads} -f {Roary_dir} -e -r -v {input.prokka_gff} &> {log}", #-e, core genes alignment using PRANK, -r, Rplots ,
        # sometimes you may need to set -i (minimum percentage identity for blastp) and -s (dont split paralogs), according to the organism
        "python /home/software/Roary/contrib/roary_plots/roary_plots.py {Roary_dir}accessory_binary_genes.fa.newick {Roary_dir}gene_presence_absence.csv && \
        mv -t {Roary_dir} pangenome_frequency.png pangenome_matrix.png pangenome_pie.png",
        """printf "\ngenerating a tree using Fasttree \n" &>> {log} && printf "\ngenerating a tree using Fasttree\n" """,
        "FastTree -nt -gtr {Roary_dir}core_gene_alignment.aln > {Roary_dir}core_gene_tree.nwk &>> {log}"
        ]
        for c in commands:
            shell(c)
