for file in *.vcf; do
	java -jar /media/perimeter/r2/srcs/snpEff_latest_core/snpEff/snpEff.jar -v eeu $file -c /media/perimeter/r2/srcs/snpEff_latest_core/snpEff/snpEff.config -no-intergenic -lof -no-downstream -no-upstream -no-intron -no-utr -nodownload -csvStats ${file}.snp.csv > tmp.snpeff.vcf
done

for file in *.vcf.snp.csv;do grep -e missense_variant -e synonymous_variant ${file} | awk '{print "'$file'""\t"$0}' >>eeu.snpeff.tab ;done
