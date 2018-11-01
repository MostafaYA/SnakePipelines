#!/bin/bash
#Mostafa.Abdel-Glil@fli.de, date: October, 24, 2018
# A bash script to find the broken links in a folder and once found, excute a copy command instead of linking
#things done
#list the broken links in

#***************************************************************NOt complete*****************************************************
#TO-DO
# summarize the broken links in a file, track the source of the original linked data and excute the cp command

directory=$1
#list all symlinks in a directory
#ls -l $directory | grep ^l
find $directory -type l -ls


for testing in $(find $directory -type l -exec sh -c 'file -b "$1" | grep -q ^broken' sh {} \; -print | xargs -n 1 basename 2> /dev/null); #find and print broken links, only the base name
#source: https://unix.stackexchange.com/questions/34248/how-can-i-find-broken-symlinks
  do  testingNr=$(echo $testing | awk '{print NF}');
  echo $testing;
  echo $testingNr;
    if [ $testingNr != 0 ];
      then cp -r /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/results/5_Prokka/16S0529/16S0529.gff ./123;
    fi ;
  done

#example output
#ERR152322.gff

#for testing in $(find . -type l -exec sh -c 'file -b "$1" | grep -q ^broken' sh {} \; -print | xargs -n 1 basename 2> /dev/null); #find and print broken links, only the base name
#source: https://unix.stackexchange.com/questions/34248/how-can-i-find-broken-symlinks
#  do  testingNr=$(echo $testing | awk '{print NF}');
#  echo $testing;
#  echo $testingNr;
#    if [ $testingNr != 0 ];
#      then cp -r /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project/results/5_Prokka/16S0529/16S0529.gff ./123;
#    fi ;
#  done
