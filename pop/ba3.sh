
#setp1: Exact Non-missing

#step2: bcf2ba3 

#step3 poptrans

#Convert to real ba3, not sure if this step is needed.

cat /media/perimeter/r2/eeu/geneflow/eeu.pop.n10.ba3 | awk '{if ($4=="A") print $1"\t"$2"\t"$3"\t"100"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="A") print $1"\t"$2"\t"$3"\t"$4"\t"100;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($4=="T") print $1"\t"$2"\t"$3"\t"110"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="T") print $1"\t"$2"\t"$3"\t"$4"\t"110;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($4=="G") print $1"\t"$2"\t"$3"\t"120"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="G") print $1"\t"$2"\t"$3"\t"$4"\t"120;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($4=="C") print $1"\t"$2"\t"$3"\t"130"\t"$5;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' | awk '{if ($5=="C") print $1"\t"$2"\t"$3"\t"$4"\t"130;else print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > eeu.pop.n10f.ba3 
