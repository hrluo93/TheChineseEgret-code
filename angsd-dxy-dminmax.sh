#!/bin/bash

#conda activate angsd

cd /media/luoub/r2/eeu/geneflow2/dxy
mkdir res
spe1=cyz
spe2=xrt
genome=/media/luoub/r2/eeu/eeu.final2.fasta

for i in chr1 chr2 chr3 chr4a chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29
do     
angsd -nThreads 8 -anc $genome -ref $genome  -r $i -bam allbam.list -skipTriallelic 1 -setMinDepth 5 -setMaxDepth 1000 -GL 1 -out $i-all -doMaf 1 -SNP_pval 1e-6 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doSaf 1 -doCounts 1 

realSFS print $i-all.saf.idx | cut -f 1-2 > $i-sites.txt
angsd sites index $i-sites.txt
angsd -nThreads 8 -anc $genome  -ref $genome -r $i -bam ${spe1}bam.list -skipTriallelic 1 -setMinDepth 5 -setMaxDepth 1000 -GL 1 -out $i-pop${spe1} -doMaf 1 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doSaf 1 -doCounts 1 -sites $i-sites.txt

angsd -nThreads 8 -anc $genome  -ref $genome -r $i -bam ${spe2}bam.list -skipTriallelic 1 -setMinDepth 5 -setMaxDepth 1000 -GL 1 -out $i-pop${spe2} -doMaf 1 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doSaf 1 -doCounts 1 -sites $i-sites.txt
angsd -nThreads 8 -anc $genome  -ref $genome -r $i -bam ${spe1}bam.list -skipTriallelic 1 -setMinDepth 5 -setMaxDepth 1000 -GL 1 -out $i-${spe1}pi -doMaf 1 -SNP_pval 1e-6 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doSaf 1 -doCounts 1
angsd -nThreads 8 -anc $genome  -ref $genome -r $i -bam ${spe2}bam.list -skipTriallelic 1 -setMinDepth 5 -setMaxDepth 1000 -GL 1 -out $i-${spe2}pi -doMaf 1 -SNP_pval 1e-6 -doMajorMinor 1 -minMapQ 20 -minQ 20 -doSaf 1 -doCounts 1

realSFS ${i}-${spe1}pi.saf.idx > ${i}-${spe1}pi.sfs
realSFS saf2theta ${i}-${spe1}pi.saf.idx -outname ${spe1}${i} -sfs ${i}-${spe1}pi.sfs
thetaStat print ${spe1}${i}.thetas.idx > ${spe1}${i}.thetas.table
awk '{print $1"\t"$2"\t"$2"\t"10^$4}' ${spe1}${i}.thetas.table > ${spe1}${i}.thetapi.bed
grep -w ${i} ${genome}.fai | cut -f1,2 > ${spe1}.${i}.g
bedtools complement -i ${spe1}${i}.thetapi.bed -g ${spe1}.${i}.g | awk '{print $1"\t"$2"\t"$3"\t"0}' > ${spe1}.${i}.invar.bed
cat ${spe1}.${i}.invar.bed ${spe1}${i}.thetapi.bed | bedtools sort -i - -g ${spe1}.${i}.g > ${spe1}${i}.allpi.bed
bedtools makewindows -g ${spe1}.${i}.g -w 50000 > ${i}.50k.bed
bedtools map -a ${i}.50k.bed -b ${spe1}${i}.allpi.bed -c 4 -o sum | awk '{print $1"\t"$2"\t"$3"\t"$4/50000}' > res/${spe1}.${i}.50kpi.bed
grep -w ${i} /media/luoub/r2/eeu/geneflow2/varfst2/eeu.rho1.bed > ${i}.rho.bed
bedtools map -a ${i}.rho.bed -b ${spe1}${i}.allpi.bed -c 4 -o sum | awk '{print $1"\t"$2"\t"$3"\t"$4/($3-$2+1)}' > res/${i}.${spe1}pi.rho.bed

realSFS ${i}-${spe2}pi.saf.idx > ${i}-${spe2}pi.sfs
realSFS saf2theta ${i}-${spe2}pi.saf.idx -outname ${spe2}${i} -sfs ${i}-${spe2}pi.sfs
thetaStat print ${spe2}${i}.thetas.idx > ${spe2}${i}.thetas.table
awk '{print $1"\t"$2"\t"$2"\t"10^$4}' ${spe2}${i}.thetas.table > ${spe2}${i}.thetapi.bed
grep -w ${i} ${genome}.fai | cut -f1,2 > ${spe2}.${i}.g
bedtools complement -i ${spe2}${i}.thetapi.bed -g ${spe2}.${i}.g | awk '{print $1"\t"$2"\t"$3"\t"0}' > ${spe2}.${i}.invar.bed
cat ${spe2}.${i}.invar.bed ${spe2}${i}.thetapi.bed | bedtools sort -i - -g ${spe2}.${i}.g > ${spe2}${i}.allpi.bed
bedtools map -a ${i}.50k.bed -b ${spe2}${i}.allpi.bed -c 4 -o sum | awk '{print $1"\t"$2"\t"$3"\t"$4/50000}' > res/${spe2}.${i}.50kpi.bed
bedtools map -a ${i}.rho.bed -b ${spe2}${i}.allpi.bed -c 4 -o sum | awk '{print $1"\t"$2"\t"$3"\t"$4/($3-$2+1)}' > res/${i}.${spe2}pi.rho.bed


gzip -d ${i}-pop${spe1}.mafs.gz -c > ${i}-pop${spe1}.mafs
gzip -d ${i}-pop${spe2}.mafs.gz -c > ${i}-pop${spe2}.mafs
len=`grep -w ${i} ${genome}.fai | cut -f2`
Rscript /media/luoub/r2/calcDxy.R -p ${i}-pop${spe1}.mafs -q ${i}-pop${spe2}.mafs -t ${len}
mv Dxy_persite.txt ${i}-Dxy_persite.txt
awk '{print $1"\t"$2"\t"$2"\t"$3}' ${i}-Dxy_persite.txt | grep -w ${i} > ${i}-Dxy_persite.bed
bedtools map -a ${i}.50k.bed -b ${i}-Dxy_persite.bed -c 4 -o min,max > res/${i}.50kdmindmax.bed
bedtools complement -i ${i}-Dxy_persite.bed -g ${spe1}.${i}.g | awk '{print $1"\t"$2"\t"$3"\t"0}' > ${i}.dxyinvar.bed
cat ${i}.dxyinvar.bed ${i}-Dxy_persite.bed | bedtools sort -g ${spe1}.${i}.g -i - > ${i}.angsddxy-all.bed
bedtools map -a ${i}.50k.bed -b ${i}.angsddxy-all.bed -c 4 -o sum | awk '{print $1"\t"$2"\t"$3"\t"$4/50000}' > res/${i}.50kdxy.bed
#grep -w ${i} /media/luoub/r2/eeu/geneflow2/varfst2/eeu.rho1.bed > ${i}.rho.bed
bedtools map -a ${i}.rho.bed -b ${i}.angsddxy-all.bed -c 4 -o sum | awk '{print $1"\t"$2"\t"$3"\t"$4/($3-$2+1)}' > res/${i}.rhodxy.bed
bedtools map -a ${i}.rho.bed -b ${i}-Dxy_persite.bed -c 4 -o min,max | awk '{print $1"\t"$2"\t"$3"\t"$4/($3-$2+1)"\t"$5/($3-$2+1)}' > res/${i}.rhodminmax.bed
done
