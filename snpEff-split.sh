for file in *.vcf; do
  mkdir ${file}
  for sample in `bcftools query -l $file`; do
    bcftools view -c1 -s $sample -o ${file}/${sample}.vcf $file
  done
done


