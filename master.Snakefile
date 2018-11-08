"""----------------------------
Mostafa.Abdel-Glil@fli.de, date: October, 15, 2018, Friedrich-Loeffler-Institut (https://www.fli.de/)
-------------------------------
# To-Do
# config file
- make 2 config files (samples and other settings)
#tools
- define which version of software to be used via specifying the full path
-------------------------------
"""
configfile: "config.yaml"
sample=config["samples"]
genome=config["genomes"]
db=config["DB"]["kraken"],
db_krona=config["DB"]["krona"]
scripts_dir=config["tools"]["scripts_dir"]
#configfile: "assembly.config.yaml"
#configfile: "config/samples.yaml" #make 2 config files (samples and other settings)

import os
import tempfile
#import glob

#TMP_DIR_ROOT = config['tmp_dir_root']

snakefiles = os.path.join(config["directories"]["snakemake_folder"], "Snakefiles/")
include: snakefiles + "defineFolders.Snakefile"
include: snakefiles + "denovoassembly.Snakefile"
include: snakefiles + "pangenome.Snakefile"
include: snakefiles + "virulenceAndAMR.Snakefile"


rule all: #unneeded rules could be deactivated by commenting the line (#)
          #deactivated runs will only run if they are needed for other rules
    input:
        #"""Rules for denovoassembly.Snakefile""",
        #data_preparation,
        FastQC,
        SICKLE,
        SPAdes,
        Kraken_fastq,
        Kraken_contig,
        Quast,
        denovo_assembly_report,
        #"""Rules for the pangenome.Snakefile""",
        Prokka,
        Prokka_assemblies,
        Roary,
        FastTree,
        #"""Rules for the virulenceAndAMR.Snakefile""",
        #dataset_link,
        aribaAMR,   ##search for AMR genes using ariba as a tool and CARD as AMR database
        aribaAMRsummary,
        abricateAMR,
        AMRabricatesummary,
        aribaVF,
        aribaVFsummary,
        abricateVF,
        VFabricatesummary,
