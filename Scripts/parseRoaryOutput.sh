#!/bin/bash
#Mostafa.Abdel-Glil@fli.de, date: October, 26, 2018, Friedrich-Loeffler-Institut (https://www.fli.de/)
#these are number of commands to parse the roary output files

#dependencies
#- seqkit

#TO-DO
# extract the 1st column that has the genes name and provide it to seqkit to extract the gene sequence from the pan_genome_reference.fa
# make a temporary output file or a pipe during determining the cutoff 
#Variables
roaryOutFolder=$1
percentageCutOff=$2

#check requirements
echo ""
if [[ -e $roaryOutFolder/gene_presence_absence.Rtab && -e $roaryOutFolder/pan_genome_reference.fa ]] ;
then
  echo "found files gene_presence_absence.Rtab and pan_genome_reference.fa" ;
elif [[ ! -e $roaryOutFolder/pan_genome_reference.fa  ]];
then
  echo "can not find the file: pan_genome_reference.fa";
  exit 1;
else
    [[ ! -e $roaryOutFolder/gene_presence_absence.Rtab ]];
    echo "can not find the gene_presence_absence.Rtab file";
    exit 1;
fi

#Parsing the gene_presence_absence.Rtab file
#1- scan the values within the file gene_presence_absence.Rtab, they should include only 0 and 1, otherwise report other values
binary=$(awk 'NR>=2' $roaryOutFolder/gene_presence_absence.Rtab | awk 'BEGIN{FS=OFS="\t"}{for(i=2;i<=NF;++i)print $i}' | sort | uniq | tr -d '\n')
#awk 'NR>=2' $roaryOutFolder/gene_presence_absence.Rtab | awk 'BEGIN{FS=OFS="\t"}{for(i=2;i<=NF;++i)print $i}' | sort | uniq | tr -d '\n'
if [[ $binary != 01 ]];
then
    echo "chracters other than 0 and 1 were detected in gene_presence_absence.Rtab, exit";
    awk 'NR>=2' $roaryOutFolder/gene_presence_absence.Rtab | awk 'BEGIN{FS=OFS="\t"}{for(i=2;i<=NF;++i)print $i}' | sort | uniq | tr '\n' ' '
    exit 1
fi
printf "\n------------------------------------\n"
#2- add 2 columns at the end of the file with the total number of genes present and the percentage of each gene
awk 'BEGIN{FS=OFS="\t"}{sum=0; cout=NF-1;
    for (i=1; i<=NF; i++) { sum+= $i }
    print (NR==1 ? $0 OFS "Total" OFS "Percentage": $0 OFS sum OFS sum/cout*100)}' $roaryOutFolder/gene_presence_absence.Rtab > 123test #To-Do: make a temporary output file or a pipe
#3- make a cutoff and print all genes equal or greater than the cutoff value
awk -v var="$percentageCutOff" '{Percentage=NF; if ($Percentage >= var)  { print ;}}' 123test #To-Do: make a temporary output file or a pipe


#old
#awk 'NR>=2' $roaryOutFolder/gene_presence_absence.Rtab | \
#  awk 'BEGIN{FS=OFS="\t"} {for (i=1; i<=NF-2; i++) $i = $(i+1); NF-=0; print}' |\
#  awk 'BEGIN{FS=OFS="\t"}{sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $0 OFS sum}'
#
#
#awk 'BEGIN{FS=OFS="\t"}{sum=0; cout=NF-1; for (i=1; i<=NF; i++) { sum+= $i } print $0 OFS sum OFS cout OFS sum/cout*100}' #the one

#NR-1 OFS $0

#in a line
#awk 'BEGIN{FS=OFS="\t"} {for (j=1; j<=NF-2; j++) $j = $(j+1); NF-=0; } {sum=0; for (i=1; i<=NF; i++) { sum+= $i } print $0 OFS sum}' gene_presence_absence.Rtab

#old
#awk 'NR>=2 {for (i=1; i<=NF-2; i++) $i = $(i+1); NF-=0; print}' $roaryOutFolder/gene_presence_absence.Rtab
