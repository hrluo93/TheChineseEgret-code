# TheChineseEgret-code
Custom scripts for the Chinese egret population genomics


create-slim-genome.py:  Generates a list of SLiM initializeGenomicElement(...) statements to simulate a random chromosome organization with realistic gene structure. This script is a Python recoding of SLiM_Manual.pdf (pp. 164–165), simulating g1(exon)/g2(intron) regions using log-normal distributions and g3 with a uniform distribution, arranged as g3→g1→(g2,g1)*→g3. The mean (in original scale) and sigma (in log space) for g1 and g2 can be estimated from real annotations using compute_g1_g2_g3_params.py and passed via --g1-mean/--g1-sigma --g2-mean/--g2-sigma -g3min/-g3max in this script. (Recommend)

create-slim-genome2.py: Deep-learning from an actual chromosome SLiM-formatted map (e.g., from gff_to_slim.py), using a simple fully connected deep neural network to learn and simulate the positional length and count distributions of genomic elements, and generate SLiM-formatted regions accordingly. TensorFlow 2.16.1 is required and tested on NVIDIA Ada Lovelace GPUs (Compute Capability 8.9).
For avian chromosomes, a reasonable fallback configuration (if no sufficient training data is available) can be:
-g3min 300 -g3max 80000 --g1-mean 165.63095 --g1-sigma 0.692356 --g2-mean 3569.4043 --g2-sigma 1.523129 ###For customized fallback parameters from real genomic data, you may use the companion script compute_g1_g2_g3_params.py to extract empirical log-normal estimates from BED/GFF-derived length tables.
