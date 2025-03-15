import numpy as np
import hicstraw
import os
import pandas as pd

##Code sourced from https://github.com/deeptools/HiCExplorer/issues/821#issuecomment-1316842070
hic_file = 'eeu.final2.fasta.hic'
cool_file = 'eeu.final2.50kb.cool'

data_type = 'observed' # (previous default / "main" data) or 'oe' (observed/expected)
normalization = "NONE"  # , VC, VC_SQRT, KR, SCALE, etc.
resolution = 50000

hic = hicstraw.HiCFile(hic_file)

assert resolution in hic.getResolutions(), \
    f"{resolution} is not part of the possible resolutions {','.join(hic.getResolutions())}"

chrom_sizes = pd.Series({chrom.name: chrom.length for chrom in hic.getChromosomes() if chrom.name != "All"})

# First write the chromosome sizes:
with open(hic.getGenomeID() + '.size', 'w') as fsize:
    for chrom in hic.getChromosomes():
        if chrom.name != "All":
            fsize.write(f"{chrom.name}\t{chrom.length}\n")
# Then write the counts in text file:
with open(cool_file.replace('.cool', ".txt"), 'w') as fo:
    for i in range(len(chrom_sizes)):
        for j in range(i, len(chrom_sizes)):
            chrom1 = chrom_sizes.index[i]
            chrom2 = chrom_sizes.index[j]
            result = hicstraw.straw(data_type, normalization, hic_file, chrom1, chrom2, 'BP', resolution)
            for k in range(len(result)):
                start1 = result[k].binX
                start2 = result[k].binY
                value = result[k].counts
                fo.write(f"{chrom1}\t{start1}\t{start1}\t{chrom2}\t{start2}\t{start2}\t{value}\n")

os.system(f"cooler load -f bg2 {hic.getGenomeID()}.size:{resolution} {cool_file.replace('.cool', '.txt')} {cool_file}")
