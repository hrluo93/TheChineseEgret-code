#Order contig
python /path/quarTeT-1.2.5/quartet_assemblymapper.py -r eeu.final2.fasta -q ul5.asm.hic.hap2.p_ctg.fa --keep --nofilter -p hap2.quar -t 32
ragtag.py splitasm hap2.quar.draftgenome.fasta > hap2.contig.fasta

#Mapping hic reads (require Haphic)
bwa index hap2.contig.fasta
bwa mem -5SP -t 32 hap2.contig.fasta eeu.repaired_R1.paired.fastq.gz eeu.repaired_R2.paired.fastq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o hap2.HiC.bam
/path/HapHiC/utils/filter_bam.py hap2.HiC.bam 1 --NM 3 --threads 14 | samtools view - -b -@ 14 -o hap2.hic.filtered.bam

#bam2juicer2.hic (require 3D-DNA matlock)
awk -f /path/3d-dna-201008/utils/generate-assembly-file-from-fasta.awk hap2.contig.fasta > hap2.contig.assembly
matlock bam2 juicer hap2.hic.filtered.bam hap2.links.mnd
sort -k2,2 -k6,6 hap2.links.mnd > hap2.sorted.links.mnd
bash /path/3d-dna-201008/visualize/run-assembly-visualizer.sh -p false hap2.contig.assembly hap2.sorted.links.mnd
