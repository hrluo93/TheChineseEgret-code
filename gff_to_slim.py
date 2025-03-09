#!/usr/bin/env python3
import argparse
import sys

def parse_attr(attr, key):
    """
    解析 GFF 属性字段，返回指定 key 的值。
    属性格式示例：ID=EEU000260.1;Parent=EEU000260
    """
    parts = attr.strip().split(";")
    for part in parts:
        part = part.strip()
        if part.startswith(key + "="):
            return part.split("=")[1]
    return None

def parse_gff(gff_file):
    """
    解析输入 GFF 文件，将 gene、mRNA 和 CDS 信息存入字典中。
    返回字典 genes: { gene_id: { 'chr':, 'strand':, 'start':, 'end':, 'mrna': [], 'cds': [] } }
    这里假定每个 gene 只有一个 mRNA，且 CDS 的 Parent 为 mRNA 的 ID。
    """
    genes = {}
    mrna_to_gene = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attr = fields
            start = int(start)
            end = int(end)
            if feature == "gene":
                gene_id = parse_attr(attr, "ID")
                if gene_id is None:
                    continue
                genes[gene_id] = {"chr": chrom, "strand": strand, "start": start, "end": end, "mrna": [], "cds": []}
            elif feature == "mRNA":
                mrna_id = parse_attr(attr, "ID")
                parent = parse_attr(attr, "Parent")
                if mrna_id is None or parent is None:
                    continue
                mrna_to_gene[mrna_id] = parent
                if parent in genes:
                    genes[parent]["mrna"].append(mrna_id)
            elif feature == "CDS":
                parent = parse_attr(attr, "Parent")
                if parent is None:
                    continue
                gene_id = mrna_to_gene.get(parent)
                if gene_id and gene_id in genes:
                    genes[gene_id]["cds"].append((start, end))
    return genes

def filter_overlapping_genes(genes):
    """
    过滤掉同一染色体上存在重叠的 gene（如果两个 gene 重叠，则两者均被移除），
    以确保后续输出的区段不会重叠。
    返回过滤后的 gene 字典。
    """
    # 按染色体分组 gene
    genes_by_chr = {}
    for gene_id, info in genes.items():
        chrom = info['chr']
        genes_by_chr.setdefault(chrom, []).append((gene_id, info['start'], info['end']))
    remove_set = set()
    for chrom, gene_list in genes_by_chr.items():
        gene_list_sorted = sorted(gene_list, key=lambda x: x[1])
        for i in range(len(gene_list_sorted)-1):
            gene_id, start, end = gene_list_sorted[i]
            next_gene_id, next_start, next_end = gene_list_sorted[i+1]
            # 如果后一个 gene 的起始位置小于等于前一个 gene 的结束位置，则重叠
            if next_start <= end:
                remove_set.add(gene_id)
                remove_set.add(next_gene_id)
    filtered_genes = { gene_id: info for gene_id, info in genes.items() if gene_id not in remove_set }
    return filtered_genes

def compute_gene_features(genes):
    """
    对每个 gene，根据其链方向对 CDS 进行排序，并计算相邻 CDS 间的 intron 区间。
    返回字典 gene_features: { gene_id: [ (feature_type, start, end), ... ] }，
    列表中依次为 CDS 与 intron（按转录顺序）。
    """
    gene_features = {}
    for gene_id, info in genes.items():
        cds_list = info.get("cds", [])
        if not cds_list:
            continue
        # 正链按 start 升序，负链按 start 降序（保证转录顺序为 5'->3'）
        if info["strand"] == "+":
            cds_sorted = sorted(cds_list, key=lambda x: x[0])
        else:
            cds_sorted = sorted(cds_list, key=lambda x: x[0], reverse=True)
        feats = []
        feats.append(("CDS", cds_sorted[0][0], cds_sorted[0][1]))
        for i in range(1, len(cds_sorted)):
            prev = cds_sorted[i-1]
            curr = cds_sorted[i]
            if info["strand"] == "+":
                intron_start = prev[1] + 1
                intron_end = curr[0] - 1
            else:
                intron_start = curr[1] + 1
                intron_end = prev[0] - 1
            if intron_end >= intron_start:
                feats.append(("intron", intron_start, intron_end))
            feats.append(("CDS", curr[0], curr[1]))
        gene_features[gene_id] = feats
    return gene_features

def parse_genome_lengths(genome_file):
    """
    解析染色体长度文件，返回两个字典：
      genome: { chrom: length }
      cumulative: { chrom: 累计偏移量 }
    累计偏移量按文件中出现顺序计算。
    """
    genome = {}
    cumulative = {}
    offset = 0
    with open(genome_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            chrom = parts[0]
            length = int(parts[1])
            genome[chrom] = length
            cumulative[chrom] = offset
            offset += length
    return genome, cumulative

def build_relative_segments_by_chrom(genes, gene_features, cumulative):
    """
    针对提供染色体信息，根据染色体对 gene 分组，
    在每条染色体内计算 gene 内部及 gene 间区域的相对坐标，
    然后加上染色体累积偏移量，生成全局连续坐标。
    返回列表 segments，每个元素为 (feature, global_start, global_end)。
    """
    segments = []
    genes_by_chr = {}
    for gene_id, info in genes.items():
        chrom = info["chr"]
        genes_by_chr.setdefault(chrom, []).append((gene_id, info))
    for chrom, gene_list in genes_by_chr.items():
        gene_list_sorted = sorted(gene_list, key=lambda x: x[1]["start"])
        offset = cumulative.get(chrom, 0)
        cum = 0  # 本染色体相对坐标
        first_gene = gene_list_sorted[0][1]
        if first_gene["start"] > 1:
            gap = first_gene["start"] - 1
            segments.append(("intergenic", offset + cum, offset + cum + gap - 1))
            cum += gap
        for i, (gene_id, info) in enumerate(gene_list_sorted):
            feats = gene_features.get(gene_id, [])
            if info["strand"] == "-":
                feats = list(reversed(feats))
            for feat in feats:
                ftype, fstart, fend = feat
                seg_length = abs(fend - fstart) + 1
                segments.append((ftype, offset + cum, offset + cum + seg_length - 1))
                cum += seg_length
            if i < len(gene_list_sorted) - 1:
                current_gene_end = info["end"]
                next_gene_start = gene_list_sorted[i+1][1]["start"]
                gap = next_gene_start - current_gene_end - 1
                if gap > 0:
                    segments.append(("intergenic", offset + cum, offset + cum + gap - 1))
                    cum += gap
    return segments

def build_relative_segments_single(genes, gene_features):
    """
    当未提供染色体长度文件时，假定所有 gene 来自同一染色体，
    按基因在基因组上的顺序计算连续相对坐标。
    """
    sorted_genes = sorted(genes.items(), key=lambda x: x[1]["start"])
    segments = []
    cum = 0
    first_gene = sorted_genes[0][1]
    if first_gene["start"] > 1:
        gap = first_gene["start"] - 1
        segments.append(("intergenic", cum, cum + gap - 1))
        cum += gap
    for i, (gene_id, info) in enumerate(sorted_genes):
        feats = gene_features.get(gene_id, [])
        if info["strand"] == "-":
            feats = list(reversed(feats))
        for feat in feats:
            ftype, fstart, fend = feat
            seg_length = abs(fend - fstart) + 1
            segments.append((ftype, cum, cum + seg_length - 1))
            cum += seg_length
        if i < len(sorted_genes) - 1:
            current_gene_end = info["end"]
            next_gene_start = sorted_genes[i+1][1]["start"]
            gap = next_gene_start - current_gene_end - 1
            if gap > 0:
                segments.append(("intergenic", cum, cum + gap - 1))
                cum += gap
    return segments

def output_slim(segments):
    """
    对所有区段按照全局起始坐标排序后，以 SLiM 脚本格式输出。
    格式示例：
      initializeGenomicElement(intergenic, 0, 2243);
      initializeGenomicElement(CDS, 2244, 2267);
    """
    segments_sorted = sorted(segments, key=lambda seg: seg[1])
    for seg in segments_sorted:
        ftype, start, end = seg
        print(f"initializeGenomicElement({ftype}, {start}, {end});")

def main():
    parser = argparse.ArgumentParser(
        description="将 GFF 文件转换为 SLiM 格式（添加 intron 与 intergenic 区域，负链转换为正链输出），"
                    "同时过滤掉重叠的 gene，并对输出坐标排序。"
    )
    parser.add_argument("-i", "--input", required=True, help="输入 GFF 文件")
    parser.add_argument("-g", "--genome", help="染色体长度文件，格式：chrom<TAB>length")
    args = parser.parse_args()

    genes = parse_gff(args.input)
    if not genes:
        sys.exit("未能解析到任何 gene 信息，请检查输入文件格式。")
    # 过滤掉存在重叠的 gene
    genes = filter_overlapping_genes(genes)
    if not genes:
        sys.exit("所有 gene 均存在重叠，未能得到非重叠 gene。")
        
    gene_features = compute_gene_features(genes)
    if args.genome:
        genome, cumulative = parse_genome_lengths(args.genome)
        segments = build_relative_segments_by_chrom(genes, gene_features, cumulative)
    else:
        segments = build_relative_segments_single(genes, gene_features)
    output_slim(segments)

if __name__ == "__main__":
    main()
