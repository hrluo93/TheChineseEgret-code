# TheChineseEgret-code
Custom scripts for the Chinese egret population genomics


create-slim-genome.py: Python-recode SLiM_Manual.pdf Page:164-165 generating a random chromosome organization.

create-slim-genome2.py: Deep-learning from an actual chromosome SLiM-formatted map (gff_to_slim.py), using a simple fully connected deep neural network to learn and simulate the positional length and count distributions of genomic elements, and generate SLiM-formatted regions based upon that information. Tensorflow 2.16.1 is required and tested on NVIDIA Ada Lovelace (Compute Capability 8.9). 
