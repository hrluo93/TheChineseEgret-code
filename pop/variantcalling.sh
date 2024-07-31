cd /media/perimeter/r2/eeu/bam
bwa index -p eeuf2 /media/perimeter/r2/eeu/eeu.final2.fasta
ls /media/perimeter/r2/eeute/Sample*/*_R1.fastq.gz | cut -d "_" -f 1-3 | while read loc;do
spname=`ls ${loc}_R1.fastq.gz | cut -d "-" -f 2`
bwa mem -R "@RG\tID:$spname\tLB:$spname\tPL:ILLUMINA\tSM:$spname" -t 32 eeuf2 ${loc}_R1.fastq.gz ${loc}_R2.fastq.gz | samtools sort -@ 18 -O BAM -o ${spname}.sort.bam
echo ${loc}_R1.fastq.gz ${spname} >>suc.samples
done
#samtools faidx eeu.final2.fasta
#java -jar /media/perimeter/r2/srcs/picard.jar CreateSequenceDictionary R=eeu.final2.fasta O=eeu.final2.dictls /media/perimeter/r2/eeu/bam/C*.sort.bam | cut -d "/" -f 7 | cut -d "." -f 1| while read spname;do 
cd /media/perimeter/r2/eeu/bam
samtools index /media/perimeter/r2/eeu/bam/${spname}.sort.bam
java -jar /media/perimeter/r2/srcs/picard.jar MarkDuplicates REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=/media/perimeter/r2/eeu/bam/${spname}.sort.bam OUTPUT=/media/perimeter/r2/eeu/bam/${spname}.sortmk.bam METRICS_FILE=/media/perimeter/r2/eeu/bam/${spname}.sortmk.g.bam.metric
samtools index /media/perimeter/r2/eeu/bam/${spname}.sortmk.bam
cd /media/perimeter/r2/eeu/vcf
java -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R eeu.final2.fasta -I /media/perimeter/r2/eeu/bam/${spname}.sortmk.bam -O /media/perimeter/r2/eeu/bam/${spname}.g.vcf.gz -ERC GVCF -ploidy 2
echo ${spname} >>suc.samples
done
ls /media/perimeter/r2/eeu/bam/X*.sort.bam | cut -d "/" -f 7 | cut -d "." -f 1| while read spname;do 
cd /media/perimeter/r2/eeu/bam
samtools index /media/perimeter/r2/eeu/bam/${spname}.sort.bam
java -jar /media/perimeter/r2/srcs/picard.jar MarkDuplicates REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=/media/perimeter/r2/eeu/bam/${spname}.sort.bam OUTPUT=/media/perimeter/r2/eeu/bam/${spname}.sortmk.bam METRICS_FILE=/media/perimeter/r2/eeu/bam/${spname}.sortmk.g.bam.metric
samtools index /media/perimeter/r2/eeu/bam/${spname}.sortmk.bam
cd /media/perimeter/r2/eeu/vcf
java -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R eeu.final2.fasta -I /media/perimeter/r2/eeu/bam/${spname}.sortmk.bam -O /media/perimeter/r2/eeu/bam/${spname}.g.vcf.gz -ERC GVCF -ploidy 2
echo ${spname} >>suc.samples
done
ls /media/perimeter/r2/eeu/bam/Z*.sort.bam | cut -d "/" -f 7 | cut -d "." -f 1| while read spname;do 
cd /media/perimeter/r2/eeu/bam
samtools index /media/perimeter/r2/eeu/bam/${spname}.sort.bam
java -jar /media/perimeter/r2/srcs/picard.jar MarkDuplicates REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=/media/perimeter/r2/eeu/bam/${spname}.sort.bam OUTPUT=/media/perimeter/r2/eeu/bam/${spname}.sortmk.bam METRICS_FILE=/media/perimeter/r2/eeu/bam/${spname}.sortmk.g.bam.metric
samtools index /media/perimeter/r2/eeu/bam/${spname}.sortmk.bam
cd /media/perimeter/r2/eeu/vcf
java -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R eeu.final2.fasta -I /media/perimeter/r2/eeu/bam/${spname}.sortmk.bam -O /media/perimeter/r2/eeu/bam/${spname}.g.vcf.gz -ERC GVCF -ploidy 2
echo ${spname} >>suc.samples
done
java -Xmx200g -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CombineGVCFs -R eeu.final2.fasta --variant gvcf.list -O combined.g.vcf.gz
java -Xmx200g -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar GenotypeGVCFs -R eeu.final2.fasta -V combined.g.vcf.gz -O combined.gt.vcf.gz
java -Xmx200g -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants -R eeu.final2.fasta -V combined.gt.vcf.gz --select-type-to-include INDEL -O combined.gt.Indel.vcf.gz
java -Xmx200g -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants -R eeu.final2.fasta -V combined.gt.vcf.gz --select-type-to-include SNP -O combined.gt.SNP.vcf.gz
java -Xmx200g -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar VariantFiltration -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O combined.filter.pahse1.vcf.gz -V combined.gt.SNP.vcf.gz
vcftools --gzvcf combined.filter.pahse1.vcf.gz --remove-filtered-all --recode --recode-INFO-all --out eeu.filtered.pahse1
bgzip eeu.filtered.pahse1.recode.vcf
tabix eeu.filtered.pahse1.recode.vcf.gz
bcftools view -i 'F_MISSING<0.15' -q 0.01:minor -m2 -M2 eeu.filtered.pahse1.recode.vcf.gz -Oz -o eeu.filtered.pahse2.vcf.gz
bcftools view -r chr1,chr2,chr3,chr4,chr4a,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29 eeu.filtered.pahse2.vcf.gz -Oz -o eeu.filtered.pahse2.autosome.vcf.gz
java -Xmx200g -jar /media/perimeter/r2/srcs/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar VariantFiltration -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O /media/perimeter/r2/eeu/vcf/combined.gt.filter.indel.vcf -V /media/perimeter/r2/eeu/vcf/combined.gt.Indel.vcf.gz
