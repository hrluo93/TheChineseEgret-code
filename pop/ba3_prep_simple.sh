#!/usr/bin/env bash
# BA3 prep (simple): vcftools real-time missing -> optional LD prunes -> chrom-proportional sample -> BA3
# deps: vcftools, bcftools, tabix, awk, sed, shuf
# optional deps: plink2 (if --prune plink), bcftools +prune plugin (if --prune bcftools)
set -euo pipefail

usage() {
cat <<'USAGE'
Usage:
  ba3_prep_simple.sh -i input.vcf.gz -m map.tsv -o outprefix [-n 30000] [--prune {none|plink|bcftools}] [--r2 0.2] [--win 50]

Inputs
  -i  输入VCF/BCF（建议bgzip + .tbi）
  -m  样本→群体映射表 (两列: <样本名正则> <群体名>)
  -o  输出前缀

Options
  -n         目标位点数（默认 30000）
  --prune    去连锁方式：none(默认) | plink | bcftools
  --r2       LD阈值 r^2（默认 0.2）
  --win      窗口大小（kb，默认 50；plink 使用  --indep-pairwise win step r2，其中 step=5）

Outputs
  outprefix.step1.nomiss.vcf.gz     # 实时无缺失 & 二等位SNP & 多态
  outprefix.step2.unlinked.vcf.gz   # 去连锁（若 --prune=none，则为 step1 的拷贝）
  outprefix.step3.subN.vcf.gz       # 按染色体比例抽样 N 位点
  outprefix.ba3.txt                 # BA3 输入 (sample pop locus a1 a2; ATGC→100/110/120/130; 缺失→0)

Notes
  * “实时无缺失”使用 vcftools --max-missing 1（F_MISS=0），不依赖 AN/AC TAG。
  * 多态通过 --mac 1 保证（去掉单态位点）。
USAGE
}

# defaults
N=30000
PRUNE="none"
R2=0.2
WINKB=50
STEP=5

# parse args
if [ $# -lt 1 ]; then usage; exit 1; fi
INPUT=""
MAP=""
OUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT="$2"; shift 2;;
    -m) MAP="$2"; shift 2;;
    -o) OUT="$2"; shift 2;;
    -n) N="$2"; shift 2;;
    --prune) PRUNE="$2"; shift 2;;
    --r2) R2="$2"; shift 2;;
    --win) WINKB="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

# checks
for t in vcftools bcftools tabix awk sed shuf; do
  command -v "$t" >/dev/null 2>&1 || { echo "Need $t in PATH"; exit 1; }
done
[[ -f "$INPUT" ]] || { echo "VCF not found: $INPUT"; exit 1; }
[[ -f "$MAP" ]]   || { echo "map.tsv not found: $MAP"; exit 1; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "[1/4] Real-time exact non-missing (vcftools F_MISS=0), biallelic SNPs, polymorphic ..."
# vcftools 直接做：去INDEL + 二等位 + 多态（MAC>=1）+ 实时无缺失
vcftools --gzvcf "$INPUT" \
  --remove-indels \
  --min-alleles 2 --max-alleles 2 \
  --mac 1 \
  --max-missing 1 \
  --recode --stdout \
| bgzip -c > "$OUT.step1.nomiss.vcf.gz"
tabix -f "$OUT.step1.nomiss.vcf.gz"

echo "[2/4] Optional LD pruning: $PRUNE"
case "$PRUNE" in
  none)
    cp "$OUT.step1.nomiss.vcf.gz" "$OUT.step2.unlinked.vcf.gz"
    cp "$OUT.step1.nomiss.vcf.gz.tbi" "$OUT.step2.unlinked.vcf.gz.tbi"
    ;;
  plink)
    command -v plink2 >/dev/null 2>&1 || { echo "plink2 not in PATH"; exit 1; }
    plink2 --vcf "$OUT.step1.nomiss.vcf.gz" \
           --double-id --allow-extra-chr --set-missing-var-ids @:# \
           --indep-pairwise $WINKB $STEP $R2 \
           --out "$TMP/pruned"
    plink2 --vcf "$OUT.step1.nomiss.vcf.gz" \
           --double-id --allow-extra-chr \
           --extract "$TMP/pruned.prune.in" \
           --recode vcf bgz --out "$TMP/step2"
    mv "$TMP/step2.vcf.gz" "$OUT.step2.unlinked.vcf.gz"
    tabix -f "$OUT.step2.unlinked.vcf.gz"
    ;;
  bcftools)
    # 需要 bcftools 安装了 +prune 插件
    bcftools +prune "$OUT.step1.nomiss.vcf.gz" -l "$R2" -w "${WINKB}kb" -Ov -o "$TMP/step2.vcf"
    bgzip -f "$TMP/step2.vcf"
    tabix -f "$TMP/step2.vcf.gz"
    mv "$TMP/step2.vcf.gz" "$OUT.step2.unlinked.vcf.gz"
    ;;
  *)
    echo "Invalid --prune '$PRUNE' (use none|plink|bcftools)"; exit 1;;
esac

echo "[3/4] Chromosome-proportional random sampling to N=$N"
# 统计可用位点数（基于将要被抽样的文件，也就是 step2）
bcftools view -H "$OUT.step2.unlinked.vcf.gz" | cut -f1 | awk '{c[$1]++} END{for(k in c) print k"\t"c[k]}' | sort -k1,1 > "$TMP/counts.tsv"
TOTAL=$(awk '{s+=$2} END{print s+0}' "$TMP/counts.tsv")

if [ "$TOTAL" -lt "$N" ]; then
  echo "Warning: only $TOTAL SNPs available; will sample all."
  N="$TOTAL"
fi

# 计算每条染色体的抽样数，精确校正到总数 N
awk -v N="$N" 'BEGIN{OFS="\t"}
{chr[NR]=$1; cnt[NR]=$2; total+=$2}
END{
  sum=0
  for(i=1;i<=NR;i++){
    share=N*cnt[i]/total
    base=int(share+0.5)
    if(base>cnt[i]) base=cnt[i]
    alloc[i]=base; sum+=base; frac[i]=share-base
  }
  need=N-sum
  if(need>0){
    while(need>0){
      best=-1; bi=-1
      for(i=1;i<=NR;i++) if(alloc[i]<cnt[i] && frac[i]>best){best=frac[i]; bi=i}
      if(bi<0) break
      alloc[bi]++; need--; frac[bi]=-1
    }
  } else if(need<0){
    while(need<0){
      worst=1e9; wi=-1
      for(i=1;i<=NR;i++) if(alloc[i]>0 && frac[i]<worst){worst=frac[i]; wi=i}
      if(wi<0) break
      alloc[wi]--; need++; frac[wi]=1e9
    }
  }
  for(i=1;i<=NR;i++) print chr[i], alloc[i]
}' "$TMP/counts.tsv" > "$TMP/targets.tsv"

# 逐染色体随机抽位点 (chrom\tpos)
: > "$TMP/pos.list"
while read -r CHR K; do
  if [ "$K" -gt 0 ]; then
    bcftools view -H -r "$CHR" "$OUT.step2.unlinked.vcf.gz" \
      | awk -v OFS="\t" '{print $1,$2}' | shuf -n "$K" >> "$TMP/pos.list"
  fi
done < "$TMP/targets.tsv"

sort -k1,1 -k2,2n "$TMP/pos.list" | awk 'NF' > "$TMP/pos.sorted.tsv"
bcftools view -T "$TMP/pos.sorted.tsv" -Oz -o "$OUT.step3.sub${N}.vcf.gz" "$OUT.step2.unlinked.vcf.gz"
tabix -f "$OUT.step3.sub${N}.vcf.gz"

echo "[4/4] VCF -> BA3 (sample  pop  locus  a1  a2; ATGC->100/110/120/130; missing->0)"
# 中间：sample \t locus \t a1 \t a2
bcftools query -f '[%SAMPLE\t%CHROM-%POS\t%TGT\n]' "$OUT.step3.sub${N}.vcf.gz" \
 | sed 's/[|/]/\t/g' | awk 'NF>=4{print $1"\t"$2"\t"$3"\t"$4}' > "$TMP/intermediate.tsv"

# 样本→群体（map.tsv 第一列为正则，第二列群体名）
awk 'BEGIN{OFS="\t"} FNR==NR{re[++n]=$1;pop[n]=$2;next}{
  p="UNASSIGNED";
  for(i=1;i<=n;i++){ if($1 ~ re[i]){ p=pop[i]; break } }
  print $1,p,$2,$3,$4
}' "$MAP" "$TMP/intermediate.tsv" > "$TMP/with_pop.tsv"

# ATGC->编码；缺失/./N->0
awk 'BEGIN{OFS="\t"}
function code(a){ if(a=="A")return 100; if(a=="T")return 110; if(a=="G")return 120; if(a=="C")return 130; if(a=="."||a=="0"||a=="N")return 0; return a }
{$4=code($4); $5=code($5); print $1,$2,$3,$4,$5}' "$TMP/with_pop.tsv" > "$OUT.ba3.txt"

NS=$(bcftools query -l "$OUT.step3.sub${N}.vcf.gz" | wc -l | tr -d ' ')
NL=$(bcftools view -H "$OUT.step3.sub${N}.vcf.gz" | wc -l | tr -d ' ')
echo "Done. Samples: $NS; Loci: $NL"
echo "Outputs:"
echo "  $OUT.step1.nomiss.vcf.gz"
echo "  $OUT.step2.unlinked.vcf.gz"
echo "  $OUT.step3.sub${N}.vcf.gz"
echo "  $OUT.ba3.txt"
