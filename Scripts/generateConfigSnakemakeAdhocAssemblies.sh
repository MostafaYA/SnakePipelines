#!/bin/bash
#Mostafa.Abdel-Glil@fli.de
# A bash script to
#part of the code was modified from https://github.com/kneubert/bacterial_assembly
#use: bash config_script1.sh fastqdir outputfile genusname
#e.g., bash config_script1.sh reads/ config,yaml genusname
# Example:config_script.sh fastq_folder

#***************************************************************NOt complete*****************************************************
#To-Do
#1- make a help MSG
#2- create error if the nunber of IDs less than half the number of the reads, wrong guessing of ids


#variables
# Path to Folder holding FASTQ-Files
fastqdir=$1
outputfile=$2
genusname=$3

#first scan the folder for the reads and create error if the reads don't match the required pattern
CheckSamples=$((ls ${fastqdir}/*{fastq,fastq.gz,fq,fq.gz,fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq )
if [[ $CheckSamples != 1 ]];
  then echo -e "\e[31mcan not guess the samples. Samples names must end with one of these fastq,fastq.gz,fq,fq.gz,fasta,fna,fa,fsa,fs,fnn , exit\e[39m";
  exit 1;
fi ;

#get the ID of the sample
ID_fasta=$(awk 'BEGIN{FS="."}{ print $1 }' <(ls ${fastqdir}/*{fastq,fastq.gz,fq,fq.gz,fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort)
printf "Guessing IDs for fasta files .....\nThe follwoing IDs are predicted for the Samples: \n${ID_fasta}\n\n"
printf "total number of samples =" && echo ${ID_fasta}| wc | awk '{print $2}'
#process the already assembled genomes first
#get the full path of fasta
if [[ $((ls ${fastqdir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | awk '{print NF}' | sort | uniq ) = 1 ]];
  then echo "genomes:" > ${outputfile}; #create the output file
    for genome in $(awk 'BEGIN{FS="."}{ print $1 }' <(ls ${fastqdir}/*{fasta,fna,fa,fsa,fs,fnn} 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort);
      do
        genome_FILE=$(realpath $(ls ${fastqdir}/${genome}.{fasta,fna,fa,fsa,fs,fnn}  2>/dev/null ) 2>/dev/null);
        echo "  "${genome}":" >> ${outputfile};
        echo "    fasta: "$genome_FILE >> ${outputfile};
        echo "    locustag: " "\""${genome}"\"" >> ${outputfile};
        echo "    genustag: " "\""${genusname}"\"" >> ${outputfile};
        echo "    taxonomy: " "\""${genusname}"\"" >> ${outputfile};
      done;
  fi
#create the output file
echo "samples:" >> ${outputfile}
#get the ID of the sample
ID=$(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*fastq 2> /dev/null | xargs -n 1 basename 2> /dev/null; ls ${fastqdir}/*fastq.gz 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort)
printf "Guessing IDs.....\nThe follwoing IDs are predicted for the Samples: \n${ID}\n\n"
printf "total number of samples =" && echo ${ID}| wc | awk '{print $2}'
#get the full path of the reads
#ls ${fastqdir}*fastq 2> /dev/null | xargs -n 1 basename ; ls ${fastqdir}*fastq.gz 2> /dev/null | xargs -n 1 basename #list directories content without showing the entire path
#realpath $(ls ${fastqdir}*${ID}*fastq 2> /dev/null ; ls ${fastqdir}*${ID}*fastq.gz 2> /dev/null)
for ID1 in $(awk 'BEGIN{FS="_"}{ print $1 }' <(ls ${fastqdir}/*fastq 2> /dev/null | xargs -n 1 basename 2> /dev/null; ls ${fastqdir}/*fastq.gz 2> /dev/null | xargs -n 1 basename 2> /dev/null) | uniq | sort);
  do
    FILES1=$(realpath $(ls ${fastqdir}/${ID1}*_R1_*.gz  2>/dev/null ||  ls ${fastqdir}/${ID1}*_R1_*.fastq  2>/dev/null ) 2>/dev/null)
    FILES2=$(realpath $(ls ${fastqdir}/${ID1}*_R2_*.gz  2>/dev/null ||  ls ${fastqdir}/${ID1}*_R2_*.fastq  2>/dev/null ) 2>/dev/null)
        NF1=$(echo $FILES1 | awk '{print NF}')
    if [[ $NF1 -lt 1 ]]
    then
        #if the reads dont match the format *_R1_*.gz, then check for *_1.fastq. if both formats are not there, then exit
        FILES1=$(realpath $(ls ${fastqdir}/${ID1}*_1.*.gz 2>/dev/null ||  ls ${fastqdir}/${ID1}*_1.fastq  2>/dev/null ) 2>/dev/null)
        FILES2=$(realpath $(ls ${fastqdir}/${ID1}*_2.*.gz 2>/dev/null ||  ls ${fastqdir}/${ID1}*_2.fastq  2>/dev/null ) 2>/dev/null)
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


echo "Output is written to " ${outputfile}
