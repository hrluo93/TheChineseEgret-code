verkko -d eeuverkko --local-cpus 32 --min-ont-length 49999 --hifi eeu.ccs.fq --nano /media/perimeter/r2/eeu/ul.fastq.gz --hic1 eeu.repaired_R1.paired.fastq.gz --hic2 eeu.repaired_R2.paired.fastq.gz
/media/perimeter/r2/srcs/hifiasm-0.19.8/hifiasm -o ul50k.asm -t32 --ul-cut 49999 --ul ul.fastq.gz --h1 eeu.repaired_R1.paired.fastq.gz --h2 eeu.repaired_R2.paired.fastq.gz eeu.ccs.fq
