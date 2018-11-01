#!/bin/bash
# mostafa.abdel-glil@fli.de
#based on this tutorial: https://github.com/BacterialCommunitiesAndPopulation/Wednesday18thMay/blob/master/Assembly_Tutorial.md
file=$1

if [ $# -lt 1 ]; then
	echo "USAGE: contigCoverage.sh myFile"
	echo "Output to STDOUT"
    exit 0
fi


awk '{a[$1]+=$3;++c[$1]}END{for(i in a)printf "%s\t%.1f\n", i, a[i]/c[i]}' $file
