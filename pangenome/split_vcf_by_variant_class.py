#!/usr/bin/env python3
"""Split a sequence-resolved VCF without changing its header or records.

Classification follows the House Finch allele-length logic at record level:
SNP has equal one-base alleles; MNP has equal 2..(threshold-1)-base alleles;
INDEL has a non-zero allele-length range below the threshold; and SV has an
allele-length range at or above the threshold. Equal-length alleles are retained
as complex SVs, matching the original rules. A multiallelic record is never
rewritten or genotype-remapped and is therefore assigned as one whole record.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from contextlib import ExitStack
from pathlib import Path
from typing import BinaryIO, Dict, Optional


CLASSES = ("SNP", "MNP", "INDEL", "SV")
IUPAC_SEQUENCE = re.compile(br"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+\Z")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split a sequence-resolved VCF into SNP, MNP, INDEL, and SV VCFs. "
            "All original header and record lines are copied unchanged."
        ),
        epilog=(
            "Multiallelic records remain intact and are classified using the "
            "length range across REF and every ALT allele."
        ),
    )
    parser.add_argument("input", type=Path, help="Input .vcf or .vcf.gz")
    parser.add_argument(
        "--output-prefix",
        type=Path,
        help="Output prefix (default: input name without .vcf[.gz])",
    )
    parser.add_argument(
        "--sv-threshold",
        type=int,
        default=50,
        help="Minimum allele-length difference assigned to SV (default: 50)",
    )
    parser.add_argument(
        "--plain-output",
        action="store_true",
        help="Write .vcf instead of the default compressed .vcf.gz",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing output files",
    )
    args = parser.parse_args()
    if args.sv_threshold < 2:
        parser.error("--sv-threshold must be at least 2")
    return args


def default_prefix(input_path: Path) -> Path:
    name = input_path.name
    if name.endswith(".vcf.gz"):
        name = name[:-7]
    elif name.endswith(".vcf"):
        name = name[:-4]
    return input_path.with_name(name)


def open_input(path: Path) -> BinaryIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rb")
    return path.open("rb")


def open_output(path: Path, compressed: bool) -> BinaryIO:
    if compressed:
        return gzip.open(path, "wb", compresslevel=6)
    return path.open("wb")


def parse_info_number(info: bytes, key: bytes) -> Optional[int]:
    prefix = key + b"="
    for item in info.split(b";"):
        if not item.startswith(prefix):
            continue
        values = []
        for value in item[len(prefix) :].split(b","):
            try:
                values.append(abs(int(value)))
            except ValueError:
                return None
        return max(values) if values else None
    return None


def classify_sequence_alleles(ref: bytes, alt: bytes, threshold: int) -> Optional[str]:
    if not IUPAC_SEQUENCE.fullmatch(ref) or alt in (b"", b"."):
        return None

    alternate_alleles = alt.split(b",")
    if not all(IUPAC_SEQUENCE.fullmatch(allele) for allele in alternate_alleles):
        return None

    lengths = [len(ref), *(len(allele) for allele in alternate_alleles)]
    bp_diff = max(lengths) - min(lengths)
    if bp_diff == 0:
        if len(ref) == 1:
            return "SNP"
        if len(ref) < threshold:
            return "MNP"
        return "SV"
    if bp_diff < threshold:
        return "INDEL"
    return "SV"


def classify_record(line: bytes, threshold: int) -> Optional[str]:
    fields = line.rstrip(b"\r\n").split(b"\t", 8)
    if len(fields) < 8:
        raise ValueError("record has fewer than 8 VCF columns")

    classification = classify_sequence_alleles(fields[3], fields[4], threshold)
    if classification is not None:
        return classification

    # Symbolic alleles are not expected here, but SVLEN permits a safe SV fallback.
    svlen = parse_info_number(fields[7], b"SVLEN")
    if svlen is not None and svlen >= threshold:
        return "SV"
    return None


def split_vcf(args: argparse.Namespace) -> Dict[str, int]:
    input_path = args.input.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"input VCF does not exist: {input_path}")

    prefix = (args.output_prefix or default_prefix(args.input)).resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    compressed = not args.plain_output
    suffix = ".vcf.gz" if compressed else ".vcf"
    final_paths = {name: Path(f"{prefix}.{name}{suffix}") for name in CLASSES}
    summary_path = Path(f"{prefix}.variant_class_counts.tsv")
    unclassified_path = Path(f"{prefix}.UNCLASSIFIED{suffix}")

    candidate_paths = [*final_paths.values(), summary_path, unclassified_path]
    if not args.force:
        existing = [str(path) for path in candidate_paths if path.exists()]
        if existing:
            raise FileExistsError(
                "output already exists; use --force to replace it: " + ", ".join(existing)
            )

    temporary = {
        name: path.with_name(f".{path.name}.tmp.{os.getpid()}")
        for name, path in {**final_paths, "UNCLASSIFIED": unclassified_path}.items()
    }
    counts = {name: 0 for name in (*CLASSES, "UNCLASSIFIED")}
    header_lines = 0
    record_lines = 0
    saw_column_header = False

    try:
        with ExitStack() as stack:
            source = stack.enter_context(open_input(input_path))
            outputs = {
                name: stack.enter_context(open_output(path, compressed))
                for name, path in temporary.items()
            }
            for line_number, line in enumerate(source, start=1):
                if line.startswith(b"#"):
                    header_lines += 1
                    saw_column_header |= line.startswith(b"#CHROM\t")
                    for output in outputs.values():
                        output.write(line)
                    continue

                record_lines += 1
                try:
                    classification = classify_record(line, args.sv_threshold)
                except ValueError as exc:
                    raise ValueError(f"line {line_number}: {exc}") from exc
                output_class = classification or "UNCLASSIFIED"
                outputs[output_class].write(line)
                counts[output_class] += 1

        if not saw_column_header:
            raise ValueError("input lacks a #CHROM VCF header line")
        if sum(counts.values()) != record_lines:
            raise RuntimeError("internal count mismatch")

        for name, final_path in final_paths.items():
            os.replace(temporary[name], final_path)
        if counts["UNCLASSIFIED"]:
            os.replace(temporary["UNCLASSIFIED"], unclassified_path)
        else:
            temporary["UNCLASSIFIED"].unlink()

        summary_tmp = summary_path.with_name(f".{summary_path.name}.tmp.{os.getpid()}")
        with summary_tmp.open("w", encoding="ascii", newline="") as summary:
            summary.write("class\trecords\n")
            for name in (*CLASSES, "UNCLASSIFIED", "TOTAL"):
                value = record_lines if name == "TOTAL" else counts[name]
                summary.write(f"{name}\t{value}\n")
            summary.write(f"# sv_threshold\t{args.sv_threshold}\n")
            summary.write(f"# header_lines\t{header_lines}\n")
        os.replace(summary_tmp, summary_path)
    except Exception:
        for path in temporary.values():
            path.unlink(missing_ok=True)
        raise

    counts["TOTAL"] = record_lines
    print(f"Input: {input_path}")
    for name in CLASSES:
        print(f"{name}: {counts[name]}\t{final_paths[name]}")
    if counts["UNCLASSIFIED"]:
        print(f"UNCLASSIFIED: {counts['UNCLASSIFIED']}\t{unclassified_path}")
    print(f"TOTAL: {record_lines}")
    print(f"Summary: {summary_path}")
    return counts


def main() -> int:
    args = parse_args()
    try:
        split_vcf(args)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
    
