#peradmixture
vcftools --gzvcf combined.filter.pahse1.vcf.gz --remove-filtered-all --recode --recode-INFO-all --out eeu.filtered.pahse1
bgzip eeu.filtered.pahse1.recode.vcf
tabix eeu.filtered.pahse1.recode.vcf.gz
bcftools view -r chr1,chr2,chr3,chr4,chr4a,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29 -i 'F_MISSING<0.15' -q 0.01:minor -m2 -M2 /media/perimeter/r2/eeu/vcf/eeu.filtered.pahse1.recode.vcf.gz -Oz -o eeu.filtered.pahse2.autosome0.01.vcf.gz
plink --make-bed --vcf eeu.filtered.pahse2.autosome.vcf.gz --recode --out eeauto --allow-extra-chr --chr-set 30
awk '{x+=1}{print $1"\tSNP"x"\t"$3"\t"$4}' eeauto.map > eeautonew.map
#change bim none inter chr name
plink --allow-extra-chr --chr-set 30 --threads 20 -vcf /media/perimeter/r2/eeu/vcf/eeu.filtered.pahse2.autosome.vcf.gz --pca 10 --out eeuautopca

#heterozygosity
for i in {3..61};do
spname=`sed -n '1p' autosome.gt.tab | cut -f $i`
num=`cat autosome.gt.tab | cut -f $i | awk '{if ($1=="0/1"||$1=="0|1") print}' | wc -l`
echo $spname "\t" $num >> smaple.autosomegt1.tab
done
