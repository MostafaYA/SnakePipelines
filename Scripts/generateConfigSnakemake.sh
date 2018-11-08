#!/bin/bash -e
#Mostafa.Abdel-Glil@fli.de, Friedrich-Loeffler-Institut (https://www.fli.de/)
# A bash script to generate a config file for Snakemake workflow in yaml format

#To-Do
#2- create error if the nunber of IDs less than half the number of the reads, wrong guessing of ids
#3- call gzip if pigz is not installed, find a way to keep the original files
#4- if the fastq files aand their

#dependancy
#pigz

#make a help MSG and pass arguments; modified from scripts available here http://stothard.afns.ualberta.ca/downloads/CCT/installation.html
PROGNAME=`basename $0`
function usage {
    echo "
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
   bash ./Scripts/generateConfigSnakemake.sh -d ./fastqReads/ -o assembly.config.yaml -g Campylobacter
"
}
function error_exit {
        echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
        exit 1
}
function remove_trailing_slash {
    string="$1"
    new_string=`echo "$string" | perl -nl -e 's/\/+$//;' -e 'print $_'`
    echo $new_string
}
while [ "$1" != "" ]; do
    case $1 in
        -d | --directory )      shift
                                directory=$1
                                ;;
        -o | --output )      shift
                                outputfile=$1
                                ;;
        -g | --genus )      shift
                                genusname=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done
if [ -z $directory ]; then
    error_exit "Please use '-d' to specify the directory path where fastq/fasta files are stored. Use '-h' for help."
fi
if [ -z $outputfile ]; then
    error_exit "Please use '-o' to specify an output config file. Use '-h' for help."
fi
if [ -z $genusname ]; then
    error_exit "Please use '-g' to specify the genus name. Use '-h' for help."
fi
#variables
fastqdir=`remove_trailing_slash "$directory"`
#fastqdir=$1
#outputfile=$2
#genusname=$3

#first scan the folder for the reads and create error if the reads don't match the required pattern
CheckSamples=$((ls ${fastqdir}/*{fastq,fastq.gz,fq,fq.gz,fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq )
if [[ $CheckSamples != 1 ]];
  then echo -e "\e[31mcan not guess the samples. Samples names must end with one of these fastq,fastq.gz,fq,fq.gz - exit\e[39m";
  exit 1;
fi ;
#create the output file
echo "samples:" > ${outputfile}
#get the ID of the sample
printf "Guessing IDs.....\n"
ID=$(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*{fastq,fq} 2> /dev/null | xargs -n 1 basename 2> /dev/null; ls ${fastqdir}/*{fastq,fq}.gz 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort)
ID_fasta=$(awk 'BEGIN{FS="."}{ print $1 }' <(ls ${fastqdir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort)
printf "The follwoing IDs are predicted for the Samples: \n${ID}\n${ID_fasta}\n\n"
printf "total number of samples to be assembled =" && echo ${ID}| wc | awk '{print $2}'
printf "total number of assemblies =" && echo ${ID_fasta}| wc | awk '{print $2}'
#Compress the fastqfiles if uncompressed
Checkonlyfastq=$((ls ${fastqdir}/*fastq 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq )
Checkonlyfq=$((ls ${fastqdir}/*fq 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq )
if [[ $Checkonlyfastq = 1 ]];
  then echo -e "found uncompressed fastq file. Creating .gz files. Original files will not be affected";
  pigz --keep ${fastqdir}/*fastq
  if [[ $Checkonlyfq = 1 ]];
    then echo -e "found uncompressed fq file. Creating .gz files. Original files will not be affected";
    pigz --keep ${fastqdir}/*fq
  fi;
fi ;
#get the full path of the reads
#ls ${fastqdir}*fastq 2> /dev/null | xargs -n 1 basename ; ls ${fastqdir}*fastq.gz 2> /dev/null | xargs -n 1 basename #list directories content without showing the entire path
#realpath $(ls ${fastqdir}*${ID}*fastq 2> /dev/null ; ls ${fastqdir}*${ID}*fastq.gz 2> /dev/null)
for ID1 in $(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*{fastq,fq}.gz 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort);
  do
    FILES1=$(realpath $(ls ${fastqdir}/${ID1}*_R1_*.gz  2>/dev/null ) 2>/dev/null)
    FILES2=$(realpath $(ls ${fastqdir}/${ID1}*_R2_*.gz  2>/dev/null ) 2>/dev/null)
        NF1=$(echo $FILES1 | awk '{print NF}')
    if [[ $NF1 -lt 1 ]]
    then
        #if the reads dont match the format *_R1_*.gz, then check for *_1.fastq. if both formats are not there, then exit
        FILES1=$(realpath $(ls ${fastqdir}/${ID1}*_1.*.gz 2>/dev/null ) 2>/dev/null)
        FILES2=$(realpath $(ls ${fastqdir}/${ID1}*_2.*.gz 2>/dev/null ) 2>/dev/null)
        NF1=$(echo $FILES1 | awk '{print NF}')
        if [ $NF1 -lt 1 ]
        then
            echo -e "\e[31mfile pattern must match *ID*_R1_*.fastq or *ID*_1.fastq. Files could also be zipped .gz\e[39m"
            exit 1
        fi
    fi

    echo "  "${ID1}":" >> ${outputfile}
    echo "    fw: "$FILES1 >> ${outputfile}
    echo "    rv: "$FILES2 >> ${outputfile}
    echo "    locustag: " "\""${ID1}"\"" >> ${outputfile}
    echo "    genustag: " "\""${genusname}"\"" >> ${outputfile}
    echo "    taxonomy: " "\""${genusname}"\"" >> ${outputfile}
done
#process the already assembled genomes if present
printf " \ngenomes:\n" >> ${outputfile}
#get the full path of fasta
if [[ $((ls ${fastqdir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq ) = 1 ]];
  then #printf " \n"; echo "genomes:" >> ${outputfile}; #create the output file
    for genome in $(awk 'BEGIN{FS="."}{ print $1 }' <(ls ${fastqdir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort);
      do
        genome_FILE=$(realpath $(ls ${fastqdir}/${genome}.{fasta,fna,fa,fsa,fs,fnn}  2>/dev/null ) 2>/dev/null);
        echo "  "
        echo "  "${genome}":" >> ${outputfile};
        echo "    fasta: "$genome_FILE >> ${outputfile};
        echo "    locustag: " "\""${genome}"\"" >> ${outputfile};
        echo "    genustag: " "\""${genusname}"\"" >> ${outputfile};
        echo "    taxonomy: " "\""${genusname}"\"" >> ${outputfile};
      done;
  fi

#Databases
echo "DB:
  minikraken: /home/mostafa.abdel/dbs/miniKraken/minikraken_20171019_8GB
  kraken: /home/DB_RAM/KrakenDB
  krona: /home/DB/Krona_Taxonomy
  CanSNPer: /home/software/CanSNPer/CanSNPerDB.db
  BLAST:
    nr:
      /home/DB/BLAST/nr/nr
    nt:
      /home/DB/BLAST/nt/nt
  ariba:
    card:
      /home/DB/Ariba/card_ariba_1_1_7
    MLST:
      /home/DB/Ariba/MLST/Arcobacter_MLST/ref_db
  MLST:
    /home/DB/ChewBBACA/Arcobacter/genes.txt" >> ${outputfile}
echo "AMR_db:  /data/AGr110/mostafa/ariba_fmt_dbs/ariba_Card
VF_db: /data/AGr110/mostafa/ariba_fmt_dbs/ariba_vfdb_full
micDATA: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/data/micData.txt

AMR_db_abricate: card
VF_db_abricate: vfdb">> ${outputfile}
echo "tools:
  scripts_dir: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/Scripts
  multiqc_bin: /home/mostafa.abdel/.local/bin
  prokka_bin: /home/software/miniconda3/bin" >> ${outputfile}
echo "directories:
  snakemake_folder: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project" >> ${outputfile}
echo "Output is written to " ${outputfile}
