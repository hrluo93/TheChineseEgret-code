cp aligned/inter_30.hic eeu.final2.fasta.hic
h1d basic dump /media/perimeter/r2/eeu/hic/eeu.final2.fasta.hic 50000 chr2 --datatype rawhic -o eeu.VC --gt /media/perimeter/r2/eeu/hic/eeu.final2.g --normalize VC --juicertool /media/perimeter/r2/eeu/hic/juicer_tools.2.20.00.jar
h1d call stripeTAD /media/perimeter/r2/eeu/hic/eeu.VC/50000/observed.VC.chr2.matrix.gz 50000 chr2 -o eeu50k.chr2.TAD --gt /media/perimeter/r2/eeu/hic/eeu.final2.g
