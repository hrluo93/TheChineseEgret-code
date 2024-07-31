##from uginx
#!/bin/bash
if [ $# = 1 ]
then
    bcftools query -f '[%SAMPLE \t %CHROM-%POS \t %TGT\n] \n' $1 | sed 's/|/\t/' | sed 's/\//\t/' | sed '/^[[:space:]]*$/d;s/[[:space:]]*$//' | awk '{$2="proxypop" OFS $2} 1' | sed '$d' | sed 's/ \./ 0/' | sed 's/ \./ 0/'
else
    printf "Usage: bcf2ba3 <vcf file>"
fi

#!/bin/bash
if [ $# = 2 ]
then
    :
else
    printf "Usage: poptrans <mapfile> <ba3 file>\n"
    exit
fi
if [ -e $1 ] && [ -e $2 ]
then
    :
elif [ ! -e $1 ]
then
    printf "$1 not found!\n"
    exit
elif [ ! -e $2 ]
then
    printf "$2 not found!\n"
fi
while IFS= read -r line
do
    if [ ! -z "$line" ]
    then
	regexp=`echo $line | awk '{ print $1 }'`
	popname=`echo $line | awk '{ print $2 }'`
#	re=^9[0-9_]+.*
	awk '/'$regexp'/ { print $1 "\t" p1 "\t" $3 "\t" $4 "\t" $5}' "p1=$popname" $2
    fi
done < "$1"

cat /media/perimeter/r2/eeu/geneflow/eeu.pop.n10.ba3 | awk '{if ($4=="A") print $1"\t"$2"\t"$3"\t"100"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="A") print $1"\t"$2"\t"$3"\t"$4"\t"100;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($4=="T") print $1"\t"$2"\t"$3"\t"110"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="T") print $1"\t"$2"\t"$3"\t"$4"\t"110;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($4=="G") print $1"\t"$2"\t"$3"\t"120"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="G") print $1"\t"$2"\t"$3"\t"$4"\t"120;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($4=="C") print $1"\t"$2"\t"$3"\t"130"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="C") print $1"\t"$2"\t"$3"\t"$4"\t"130;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > eeu.pop.n10f.ba3 
