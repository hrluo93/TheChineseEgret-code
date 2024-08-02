cp aligned/inter_30.hic eeu.final2.fasta.hic
h1d basic dump /media/perimeter/r2/eeu/hic/eeu.final2.fasta.hic 50000 chr2 --datatype rawhic -o eeu.VC --gt /media/perimeter/r2/eeu/hic/eeu.final2.g --normalize VC --juicertool /media/perimeter/r2/eeu/hic/juicer_tools.2.20.00.jar
h1d call stripeTAD /media/perimeter/r2/eeu/hic/eeu.VC/50000/observed.VC.chr2.matrix.gz 50000 chr2 -o eeu50k.chr2.TAD --gt /media/perimeter/r2/eeu/hic/eeu.final2.g

#manually covert to perform the robust test in HiCExplorer
python /media/perimeter/r2/eeu/hic/hic2mcool.py
hicConvertFormat --chromosomeSizes /media/perimeter/r2/eeu/hic/eeu.final2.g --matrices /media/perimeter/r2/eeu/hic/eeu.final2.50kb.txt --inputFormat 2D-text --outputFormat h5 --outFileName eeu.final2.50k.h5 --resolutions 50000
hicCorrectMatrix diagnostic_plot --matrix /media/perimeter/r2/eeu/hic/eeu.final2.50k.h5 -o h5_diag.png
hicCorrectMatrix correct -m /media/perimeter/r2/eeu/hic/eeu.final2.50k.h5 --filterThreshold -3 5 -o eeu.final2.50k.corrected.h5
hicFindTADs -m eeu.final2.50k.corrected.h5 --thresholdComparisons 0.05 --delta 0.01 --outPrefix correctedchr2.tad --numberOfProcessors 16 --correctForMultipleTesting fdr --chromosomes chr2
hicDetectLoops -m eeu.final2.50k.corrected.h5 --pValuePreselection 0.05 --pValue 0.05 -o chr2.loop.bed --chromosomes chr2




