# TheChineseEgret-code
Custom scripts for the Chinese egret population genomics


create-slim-genome.py:  Generates a list of SLiM initializeGenomicElement(...) statements to simulate a random chromosome organization with realistic gene structure. This script is a Python recoding of SLiM_Manual.pdf (pp. 164–165), simulating g1(exon)/g2(intron) regions using log-normal distributions and g3 with a uniform distribution, arranged as g3→g1→(g2,g1)*→g3. The mean (in original scale) and sigma (in log space) for g1 and g2 can be estimated from real annotations using compute_g1_g2_g3_params.py and passed via --g1-mean/--g1-sigma --g2-mean/--g2-sigma -g3min/-g3max in this script.

create-slim-genome2.py: Deep-learning from an actual chromosome SLiM-formatted map (gff_to_slim.py), using a simple fully connected deep neural network to learn and simulate the positional length and count distributions of genomic elements, and generate SLiM-formatted regions based upon that information. Tensorflow 2.16.1 is required and tested on NVIDIA Ada Lovelace (Compute Capability 8.9). 
