vcftools --gzvcf combined.filter.pahse1.vcf.gz --remove-filtered-all --recode --recode-INFO-all --out eeu.filtered.pahse1
bgzip eeu.filtered.pahse1.recode.vcf
tabix eeu.filtered.pahse1.recode.vcf.gz
vcftools --gzvcf eeu.filtered.pahse1.recode.vcf.gz --extract-FORMAT-info GT --out all.GT.tab
bcftools view -r chr1,chr2,chr3,chr4,chr4a,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29 -i 'F_MISSING<0.15' -q 0.05:minor -m2 -M2 /media/perimeter/r2/eeu/vcf/eeu.filtered.pahse1.recode.vcf.gz -Oz -o eeu.filtered.pahse2.autosome0.05.vcf.gz
plink --make-bed --vcf eeu.filtered.pahse2.autosome0.05.vcf.gz --recode --out eeauto005 --allow-extra-chr --chr-set 30
awk '{x+=1}{print $1"\tSNP"x"\t"$3"\t"$4}' eeauto005.map > eeauto055-1.map
#change bim none inter chr name

  
