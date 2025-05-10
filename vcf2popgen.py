#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
示例脚本: 从VCF读取并输出 BayeScan 和 Genepop 格式。

使用方法:
  python vcf2popgen.py -i /path/to/input.vcf \
                               -o prefix \
                               -pop /path/to/pop.csv

结果:
  prefix.bayescan
  prefix.popgen
"""

import argparse
import vcf2popgen

def main():
    parser = argparse.ArgumentParser(
        description="Convert VCF to BayeScan and Genepop using vcf2popgen"
    )

    parser.add_argument('-i', '--input', required=True,
                        help="Input VCF file path")
    parser.add_argument('-o', '--output', required=True,
                        help="Output prefix for the resulting files (e.g. 'cyz')")
    parser.add_argument('-pop', '--popfile', required=True,
                        help="Population CSV file (e.g. 'cyzpop.csv')")

    args = parser.parse_args()

    # 读取 VCF 并关联群体分配文件
    data = vcf2popgen.read(args.input, args.popfile)

    # 输出文件，以 -o 指定的前缀 (prefix) 命名
    data.to_bayescan(f"{args.output}.bayescan")
    data.to_genepop(f"{args.output}.popgen")

    print(f"[INFO] Done! Generated files:\n  {args.output}.bayescan\n  {args.output}.popgen")

if __name__ == '__main__':
    main()
