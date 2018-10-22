#!/bin/bash
#Mostafa.Abdel-Glil@fli.de
# A bash script to
#part of the code was modified from https://github.com/kneubert/bacterial_assembly
#use: bash config_script1.sh fastqdir outputfile genusname
#e.g., bash config_script1.sh reads/ config,yaml genusname
# Example:config_script.sh fastq_folder

#To-Do
#1- make a help MSG
#2- first scan the folder for the reads and create error if the reads don't match the required pattern
#3- create error if the nunber of IDs less than half the number of the reads


#variables
# Path to Folder holding FASTQ-Files
fastqdir=$1
outputfile=$2
genusname=$3

echo "samples:" > ${outputfile}
#get the ID of the sample
ID=$(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*fastq 2> /dev/null | xargs -n 1 basename; ls ${fastqdir}/*fastq.gz 2> /dev/null | xargs -n 1 basename) | uniq | sort)
printf "Guessing IDs.....\nThe follwoing IDs are predicted for the Samples: \n${ID}\n\n"
printf "total number of samples =" && echo ${ID}| wc | awk '{print $2}'
#get the full path of the reads
#ls ${fastqdir}*fastq 2> /dev/null | xargs -n 1 basename ; ls ${fastqdir}*fastq.gz 2> /dev/null | xargs -n 1 basename #list directories content without showing the entire path
#realpath $(ls ${fastqdir}*${ID}*fastq 2> /dev/null ; ls ${fastqdir}*${ID}*fastq.gz 2> /dev/null)
for ID1 in $(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*fastq 2> /dev/null | xargs -n 1 basename; ls ${fastqdir}/*fastq.gz 2> /dev/null | xargs -n 1 basename) | uniq | sort);
  do
    FILES1=$(realpath $(ls ${fastqdir}/${ID1}*_R1_*.gz  2>/dev/null ||  ls ${fastqdir}/${ID1}*_R1_*.fastq  2>/dev/null ) 2>/dev/null)
    FILES2=$(realpath $(ls ${fastqdir}/${ID1}*_R2_*.gz  2>/dev/null ||  ls ${fastqdir}/${ID1}*_R2_*.fastq  2>/dev/null ) 2>/dev/null)
        NF1=$(echo $FILES1 | awk '{print NF}')
    if [[ $NF1 -lt 1 ]]
    then
        FILES1=$(realpath $(ls ${fastqdir}/${ID1}*_1.*.gz 2>/dev/null ||  ls ${fastqdir}/${ID1}*_1.fastq  2>/dev/null ) 2>/dev/null)
        FILES2=$(realpath $(ls ${fastqdir}/${ID1}*_2.*.gz 2>/dev/null ||  ls ${fastqdir}/${ID1}*_2.fastq  2>/dev/null ) 2>/dev/null)
        NF1=$(echo $FILES1 | awk '{print NF}')
        if [ $NF1 -lt 1 ]
        then
            echo "file pattern must match *ID*_R1_*.fastq or *ID*_1.fastq. Files could also be zipped .gz"
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

#Databases
#echo $FILES1
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
  krona:
    /home/DB/Krona_Taxonomy
  ariba:
    card:
      /home/DB/Ariba/card_ariba_1_1_7
    MLST:
      /home/DB/Ariba/MLST/Arcobacter_MLST/ref_db
  MLST:
    /home/DB/ChewBBACA/Arcobacter/genes.txt" >> ${outputfile}

echo "tools:
  scripts_dir: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/Scripts
  multiqc_bin: /home/mostafa.abdel/.local/bin" >> ${outputfile}

echo "Output is written to " ${outputfile}
