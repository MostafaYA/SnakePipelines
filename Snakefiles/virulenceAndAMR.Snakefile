"""----------------------------
Mostafa.Abdel-Glil@fli.de, date: October, 15, 2018, Friedrich-Loeffler-Institut (https://www.fli.de/)
-------------------------------
The pipeline includes the following
- ARIBA AMR screnning  - screen fastq reads for Antimicrobial resistance genes including variants (e.g, mutation) between the reads and the reference sequences
- ARIBA VF screnning - screen fastq reads for virulence factors genes including variants
- ABRicate AMR screnning - screen the fasta files for antimicrobial resistance genes. Mutation will not be reported
- ABRicate VF screnning - screen the fasta files for virulence factors genes.
- Summary of all results
-------------------------------
DONE:
- AMR screening using ariba and desired database #fastq
- summary AMR from ariba, used for visualization in phandango
- VF screening using ariba and desired database #fastq
- summary VF from ariba
- VF screening using ABRicate and desired database #assemblies

# To-Do
#addrules
- plot reulsts with the MIC data # micplot, you need to correct the name in the summary files
        - correct the names in the summary files
- plot all AMR results in R
- plot all VF results in R
-------------------------------
"""
#Rules output files
#dataset_link=expand(link_dir + "{sample}_1.fastq.gz", sample=sample),
aribaAMR=expand(ariba_amr_dir + "{sample}/report.tsv", sample=sample),
aribaAMRsummary= ariba_amr_dir + "amr_summary.csv"
abricateAMR=expand(abricate_amr_dir + "{genome}_abricate_report.tsv", genome=genome)
AMRabricatesummary= abricate_amr_dir + "amr_abricate_summary.csv"
aribaVF=expand(ariba_vf_dir + "{sample}/report.tsv", sample=sample),
aribaVFsummary= ariba_vf_dir + "vf_summary.csv"
abricateVF=expand(abricate_vf_dir + "{genome}_abricateVF_report.tsv", genome=genome)
VFabricatesummary= abricate_vf_dir + "vf_abricate_summary.csv",

"""
Raw folder contin a symbolic links for the raw reads

rule dataset_link: #this is to simplify the names of the raw data, so all samples  will be named as sample ID + _1 or _2 for fw and rv reads respect.
    input:
        fw= lambda wildcards: config["samples"][wildcards.sample]["fw"], #fw= lambda wildcards: expand(config["samples"][wildcards.sample]["fw"], sample=sample),
        rv= lambda wildcards: config["samples"][wildcards.sample]["rv"],
    output:
        raw_1=temp(link_dir + "{sample}_1.fastq.gz"),
        raw_2=temp(link_dir + "{sample}_2.fastq.gz")
    run:
        commands = [
            """"""if  ln -s {input.fw} {output.raw_1} 2> /dev/null && ln -s {input.rv} {output.raw_2} 2> /dev/null ; then \
                echo "sympolic links created"; \
                else echo "sympolic links are not supported. Copying the data...." & cp -r  {input.fw} {output.raw_1} && cp -r {input.rv} {output.raw_2} ; fi""""""
            #if the creation of symbolic links is not supported then copy the files
        ]
        for c in commands:
            shell(c)"""
"""
AMR using Ariba and AMR db
"""
rule AMRariba: #search AMR genes using ariba and an AMR database
                   #search using the fastq files
    input:
        AMR_db=config['AMR_db'],
        r1 = link_dir + "{sample}_1.fastq.gz",
        r2 = link_dir + "{sample}_2.fastq.gz"
    output:
        ariba_report = ariba_amr_dir + "{sample}/report.tsv",
    threads: 32
    benchmark: benchmarks_folder + "{sample}/ariba.txt"
    log: log_folder + "{sample}/ariba.log"
    conda: envs_folder + "ariba.yaml" #needs revision
    params:
        #ariba_env = config['ariba_env'],
        outdir = ariba_amr_dir + "{sample}"
    shell:
        "ariba run --verbose --force {input.AMR_db} {input.r1} {input.r2} {params.outdir} &> {log}"
"""
Summary AMR using Ariba and AMR db
"""
rule AMRsummary:
    input:
        ariba_report = expand(ariba_amr_dir + "{sample}/report.tsv", sample=sample)
    output:
        ariba_summary = ariba_amr_dir + "amr_summary.csv",
    conda: envs_folder + "ariba.yaml" #needs revision
    params:
        outfile_name = ariba_amr_dir + "amr_summary",
    shell:
        "ariba summary {params.outfile_name} {input.ariba_report}"
"""
AMR using ABRicate and AMR db
"""
rule AMRabricate: #search AMR genes using abricate and an AMR database
                   #search using the fasta files
    input:
        adhoc_assemblies=lambda wildcards: config["genomes"][wildcards.genome]["fasta"]
    output:
        abricate_report = abricate_amr_dir + "{genome}_abricate_report.tsv",
    threads: 32
    benchmark: benchmarks_folder + "{genome}/abricate.txt"
    #log: log_folder + "{sample}/abricate.log"
    conda: envs_folder + "abricate.yaml" #needs revision
    params:
        #ariba_env = config['ariba_env'],
        AMR_abricate=config['AMR_db_abricate'],
        #outdir = abricate_amr_dir + "{genome}"
    shell:
        "abricate --nopath --threads {threads} --db {params.AMR_abricate} {input.adhoc_assemblies} >  {output.abricate_report}"
"""
Summary AMR using Abricate and AMR db
"""
rule AMRabricatesummary:
    input:
        abricate_report = expand(abricate_amr_dir + "{genome}_abricate_report.tsv", genome=genome)
    output:
        abricate_summary = abricate_amr_dir + "amr_abricate_summary.csv",
    conda: envs_folder + "abricate.yaml" #needs revision
    params:
        #ariba_env = config['ariba_env'],
    shell:
        "abricate --summary {input.abricate_report} >> {output.abricate_summary} "
"""
Virulence factors screening using Ariba and VF db
"""
rule VFariba: #search AMR genes using ariba and an AMR database
                   #search using the fastq files
    input:
        VF_db=config['VF_db'],
        r1 = link_dir + "{sample}_1.fastq.gz",
        r2 = link_dir + "{sample}_2.fastq.gz"
    output:
        ariba_report = ariba_vf_dir + "{sample}/report.tsv",
    threads: 32
    benchmark: benchmarks_folder + "{sample}/ariba_vf.txt"
    log: log_folder + "{sample}/ariba_vf.log"
    conda: envs_folder + "ariba.yaml" #needs revision
    params:
        #ariba_env = config['ariba_env'],
        outdir = ariba_vf_dir + "{sample}"
    shell:
        "ariba run --verbose --force {input.VF_db} {input.r1} {input.r2} {params.outdir} &> {log}"
"""
Summary Virulence factors screening using Ariba and AMR db
"""
rule VFaribaSummary:
    input:
        ariba_report = expand(ariba_vf_dir + "{sample}/report.tsv", sample=sample)
    output:
        ariba_summary = ariba_vf_dir + "vf_summary.csv",
    conda: envs_folder + "ariba.yaml" #needs revision
    params:
        outfile_name = ariba_vf_dir + "vf_summary",
    shell:
        "ariba summary {params.outfile_name} {input.ariba_report}"
"""
AMR using ABRicate and AMR db
"""
rule VFabricate: #search AMR genes using abricate and an AMR database
                   #search using the fasta files
    input:
        adhoc_assemblies=lambda wildcards: config["genomes"][wildcards.genome]["fasta"]
    output:
        abricate_report = abricate_vf_dir + "{genome}_abricateVF_report.tsv",
    threads: 32
    benchmark: benchmarks_folder + "{genome}/abricateVF.txt"
    #log: log_folder + "{sample}/abricate.log"
    conda: envs_folder + "abricate.yaml" #needs revision
    params:
        #ariba_env = config['ariba_env'],
        VF_abricate=config['VF_db_abricate'],
        #outdir = abricate_amr_dir + "{genome}"
    shell:
        "abricate --nopath --threads {threads} --db {params.VF_abricate} {input.adhoc_assemblies} >  {output.abricate_report}"
"""
Summary AMR using Abricate and AMR db
"""
rule VFabricatesummary:
    input:
        abricate_report = expand(abricate_vf_dir + "{genome}_abricateVF_report.tsv", genome=genome)
    output:
        abricate_summary = abricate_vf_dir + "vf_abricate_summary.csv",
    conda: envs_folder + "abricate.yaml" #needs revision
    params:
        #ariba_env = config['ariba_env'],
    shell:
        "abricate --summary {input.abricate_report} >> {output.abricate_summary} "


"""
#commnads
source activate ariba
mkdir -p ariba_dbs
ariba getref card ariba_dbs/out.card
source deactivate
ariba prepareref -f ariba_dbs/out.card.fa -m ariba_dbs/out.card.tsv ariba_dbs/card_aribaprepareref
ariba run ariba_dbs/ariba_Card results/0_fastq_gz/SRR4451713_1.fastq.gz results/0_fastq_gz/SRR4451713_2.fastq.gz results/aribaCard

ariba getref srst2_argannot ariba_dbs/srst2_argannot

ariba getref vfdb_full ariba_dbs/vfdb_full
ariba prepareref -f ariba_dbs/vfdb_full.fa -m ariba_dbs/vfdb_full.tsv ariba_dbs/ariba_vfdb_full
ariba prepareref -f ariba_dbs/vfdb_full.fa -m ariba_dbs/vfdb_full.tsv ariba_dbs/ariba_vfdb_full_1 --threads 32 --verbose

##
        #"set +u; {params.ariba_env}; set -u &\
        #mkdir -p ariba_dbs &\
        #ariba getref card ariba_dbs/out.card &\
        #source deactivate ; ariba prepareref --force -f ariba_dbs/out.card.fa -m ariba_dbs/out.card.tsv ariba_dbs/ariba_card --threads {threads}&\
        #" for i in $(realpath {input.AMR_db} 2> /dev/null | xargs -n 1 basename 2> /dev/null); do  ariba run --verbose --force {input.AMR_db} {input.r1} {input.r2} {params.outdir}_$i &> {log}; done"

##
ariba micplot /data/AGr110/mostafa/ariba_fmt_dbs/ariba_Card/ --interrupted Azithromycin /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/data/micData.txt \
/home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/results/9_ariba_amr_dir/amr_summary1.csv micplot.AZMknowngroups
"""
