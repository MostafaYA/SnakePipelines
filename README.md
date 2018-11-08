## Introduction
------------------------------
This is a bioinformatics workflow that built using snakemake with the aim to automatic downstream processing of Next Generation Sequencing (NGS) reads, currently supports only paired-end Illumina sequence data. This project bundles a number of snakefiles for a _de novo_ assembly, pan-genome and detecting genes for antimicrobial resistance as well as virulence. This is a brief description on the snakefiles     

#### _De novo_ assembly snakefile 
`denovoassembly.Snakefile` 
- Quality assessment of raw reads using FASTQC
- quality trimming using SICKLE 
- SPAdes for a _de novo_ assembly
- Contig filteration, including lengths and coverage  
- Taxonomic classification of fastq reads using minikraken
- Taxonomic classification of assembled contigs using kraken
- Calculate read coverage stats for mapped reads
- Quality assessment of assemblies using quast
- A final informative .html report using MULTIQC

#### Pan-genome snakefile  
`pangenome.Snakefile`
- Prokka for a rapid contig annotation
- Roary to construct a pangenome
- FastTree to create a ML phylogeny

### To-Do
* replace `run` in rules with `shell`, so conda packages will be dowanloaded and used 
* allow each of the snakefiles to run independently if needed   

### Current limitation
* Visualising the phylogentic tree and plotting metdata is done outside Snakemake  
* The option `--use-conda` within snakemake is not currently feasible 
* The config file **MUST** include both fastq and fasta files 

### Dependencies 
------------------------------
The workflow was used with the following versions of software   
* snakemake v5.3.0 https://snakemake.readthedocs.io/en/stable/
* Fastqc v0.11.5 https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 
* sickle version 1.33 https://github.com/najoshi/sickle 
* SPAdes v3.12.0 http://cab.spbu.ru/software/spades/ 
* Kraken v0.10.6 https://ccb.jhu.edu/software/kraken/
* KronaTools 2.7 https://github.com/marbl/Krona
* QualiMap v.2.2.1 http://qualimap.bioinfo.cipf.es/
* QUAST v4.3 http://bioinf.spbau.ru/quast
* multiqc v1.5 https://multiqc.info/
* Prokka v1.13.3 https://github.com/tseemann/prokka
* Roary v3.6.1 https://sanger-pathogens.github.io/Roary/
* Bowtie2 v2.3.0 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* BWA v0.7.12-r1039 http://bio-bwa.sourceforge.net/
* SAMTools v1.4-14-g90b995f  http://samtools.sourceforge.net/
* FastTree v2.1.9 http://www.microbesonline.org/fasttree/

## Getting Started
------------------------------
### Setting up a project folder and obtain the latest version of the workflow 
* Set up a project folder for the run 

        mkdir NGS-Project
        cd NGS-Project


* Download the latest version from gitlab  
        
        git clone https://gitlab.com/Mostafa.Abdel-Glil/snakepipelines_bacterialgenomes.git

### Generate a config file for the pipeline 
A bash script `generateConfigSnakemake.sh ` is written to automatically generate a config file in yaml format providing a folder that holds the raw data as well as the already assembled genomes. The produced config file list all raw data and assembled genome and contain the paths of databases and scripts. Some editing to config file is essential to set up the paths for databases 
```
USAGE:
   bash ./Scripts/generateConfigSnakemake.sh -d DIR/ -o File -g STRING

DESCRIPTION:
   Generate Config file for Snakemake workflow in yaml format.

REQUIRED ARGUMENTS:
   -d, --directory DIR
      directory path where fastq/fasta files are stored.
   -o, --output STRING
      The output config file for Snakemake workflow in yaml format.
   -g, --genus STRING
      The name of genus for the fastq/fasta files

OPTIONAL ARGUMENTS:
   -h, --help
      Show this message.

EXAMPLE:
   bash ./Scripts/generateConfigSnakemake.sh -d ./fastqReads/ -o config.yaml -g Campylobacter

```
**The following paths for databases should be adjusted in the config file.** 
```
DB:
  minikraken: /home/mostafa.abdel/dbs/miniKraken/minikraken_20171019_8GB
  kraken: /home/DB_RAM/KrakenDB
  krona: /home/DB/Krona_Taxonomy
AMR_db:  /data/AGr110/mostafa/ariba_fmt_dbs/ariba_Card
VF_db: /data/AGr110/mostafa/ariba_fmt_dbs/ariba_vfdb_full
micDATA: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/data/micData.txt
AMR_db_abricate: card
VF_db_abricate: vfdb
tools:
  scripts_dir: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/Scripts
  multiqc_bin: /home/mostafa.abdel/.local/bin
directories:
  snakemake_folder: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project
```

### Running the snake pipeline 
It is always a good idea to display what the workflow will do without execution. For doing that, we will use the follwoing command. 

    snakemake -np --quiet --snakefile master.Snakefile

Execute the commands in the pipeline by removing the `-np` option  
    
    snakemake --snakefile master.Snakefile


    
#### Contact 
___

Comments should be addressed to Mostafa.Abdel-Glil@fli.de 

