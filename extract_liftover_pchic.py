#!/usr/bin/env python3
"""Reproduce FOXC2 promoter Capture Hi-C extraction and hg19->hg38 liftover.

Input: Javierre et al. (2016) PCHiC_peak_matrix_cutoff5.txt.gz, hg19.
The source matrix has 1-based inclusive HindIII fragment coordinates and
CHiCAGO scores in cell-type columns (EP = endothelial precursors).

Output: BEDPE files with 0-based, half-open anchors and ORIGINAL CHiCAGO scores.
Requires: Python >= 3.9 and the UCSC liftOver executable + hg19ToHg38 chain.
No R/CHiCAGO rerun, RDS, or design archive is required for this workflow.
"""

import argparse
import csv
import gzip
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

REQUIRED_COLUMNS = (
    "baitChr", "baitStart", "baitEnd", "baitID", "baitName",
    "oeChr", "oeStart", "oeEnd", "oeID",
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--matrix", type=Path, required=True,
                   help="Javierre peak matrix (.txt or .txt.gz), hg19")
    p.add_argument("--chain", type=Path, required=True,
                   help="UCSC hg19ToHg38.over.chain.gz")
    p.add_argument("--liftover", default="liftOver",
                   help="Path to or command name of UCSC liftOver")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--bait-id", required=True,
                   help="Exact hg19 promoter bait ID, e.g. 278067 for FOXC2")
    p.add_argument("--bait-chrom", default="16",
                   help="Chromosome of hg19 bait (default: 16)")
    p.add_argument("--cell-column", default="EP",
                   help="CHiCAGO score column (EP = endothelial precursors)")
    p.add_argument("--min-score", type=float, default=5.0)
    p.add_argument("--prefix", default="FOXC2_EndP")
    p.add_argument("--expected-selected", type=int, default=None,
                   help="Optional sanity check for number of extracted interaction pairs")
    p.add_argument("--expected-hg38-bait", default=None,
                   help="Optional sanity check: chr16:86558346-86569633 (BED coordinates)")
    return p.parse_args()


def chrom(c):
    c = c.strip()
    if not c:
        raise ValueError("Missing chromosome name")
    return c if c.startswith("chr") else "chr" + c


def check_interval(start_str, end_str, rownum, label):
    start, end = int(start_str), int(end_str)
    if start < 1 or end < start:
        raise ValueError(f"Bad 1-based {label} interval in input line {rownum}: {start}-{end}")
    # Input rmap coordinates are 1-based inclusive. Convert to BED 0-based half-open.
    return str(start - 1), str(end)


def extract(matrix, out_bedpe, bait_bed, oe_bed, args):
    opener = gzip.open if matrix.suffix == ".gz" else open
    count = 0
    with (opener(matrix, "rt", newline="") as source,
          out_bedpe.open("w", newline="") as bedpe_f,
          bait_bed.open("w", newline="") as bait_f,
          oe_bed.open("w", newline="") as oe_f):
        reader = csv.DictReader(source, delimiter="\t")
        missing = [c for c in (*REQUIRED_COLUMNS, args.cell_column)
                   if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Missing input columns: {missing}. Header={reader.fieldnames}")
        bedpe_writer = csv.writer(bedpe_f, delimiter="\t", lineterminator="\n")
        bait_writer = csv.writer(bait_f, delimiter="\t", lineterminator="\n")
        oe_writer = csv.writer(oe_f, delimiter="\t", lineterminator="\n")
        for line_num, row in enumerate(reader, start=2):
            if row["baitID"].strip() != args.bait_id:
                continue
            if chrom(row["baitChr"]) != chrom(args.bait_chrom):
                continue
            raw_score = row[args.cell_column].strip()
            if raw_score in ("", "NA", "NaN", "nan"):
                continue
            score = float(raw_score)
            if score < args.min_score:
                continue
            count += 1
            uid = f"{args.prefix}_{count}"
            bstart, bend = check_interval(row["baitStart"], row["baitEnd"], line_num, "bait")
            ostart, oend = check_interval(row["oeStart"], row["oeEnd"], line_num, "other end")
            b = [chrom(row["baitChr"]), bstart, bend]
            o = [chrom(row["oeChr"]), ostart, oend]
            bedpe_writer.writerow(b + o + [uid, raw_score, ".", "."])
            bait_writer.writerow(b + [uid])
            oe_writer.writerow(o + [uid])
    if count == 0:
        raise ValueError("No interactions selected; check bait ID, chromosome, score column/threshold")
    if args.expected_selected is not None and count != args.expected_selected:
        raise ValueError(f"Selected {count}; expected {args.expected_selected}. Check input matrix/version")
    return count


def run_liftover(exe, input_bed, chain, converted, unmapped):
    subprocess.run([exe, str(input_bed), str(chain), str(converted), str(unmapped)],
                   check=True)


def load_mapped(path):
    mapped = defaultdict(list)
    with path.open(newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 4:
                raise ValueError(f"Invalid BED4 in {path}: {row}")
            chr_, start, end, uid = row[:4]
            if int(start) >= int(end):
                raise ValueError(f"Invalid mapped interval in {path}: {row}")
            mapped[uid].append([chr_, start, end])
    return mapped


def reconstruct(original, baits, oes, out_all, out_cis, failed_path):
    counts = Counter()
    hg38_baits = set()
    with (original.open(newline="") as in_f,
          out_all.open("w", newline="") as all_f,
          out_cis.open("w", newline="") as cis_f,
          failed_path.open("w", newline="") as failed_f):
        reader = csv.reader(in_f, delimiter="\t")
        all_writer = csv.writer(all_f, delimiter="\t", lineterminator="\n")
        cis_writer = csv.writer(cis_f, delimiter="\t", lineterminator="\n")
        fail_writer = csv.writer(failed_f, delimiter="\t", lineterminator="\n")
        fail_writer.writerow(["interaction_id", "bait_mappings", "other_end_mappings", "reason"])
        for row in reader:
            if not row:
                continue
            counts["selected"] += 1
            uid = row[6]
            bait_hits, oe_hits = baits.get(uid, []), oes.get(uid, [])
            if len(bait_hits) != 1 or len(oe_hits) != 1:
                counts["failed"] += 1
                fail_writer.writerow([uid, len(bait_hits), len(oe_hits),
                                      "require exactly one mapped interval for each end"])
                continue
            bait, oe = bait_hits[0], oe_hits[0]
            converted = bait + oe + row[6:]
            all_writer.writerow(converted)
            hg38_baits.add(":".join((bait[0], bait[1] + "-" + bait[2])))
            counts["converted"] += 1
            if bait[0] == oe[0]:
                cis_writer.writerow(converted)
                counts["cis"] += 1
    return counts, sorted(hg38_baits)


def main():
    args = parse_args()
    if not args.matrix.is_file():
        raise FileNotFoundError(args.matrix)
    if not args.chain.is_file():
        raise FileNotFoundError(args.chain)
    args.outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.prefix}_score{args.min_score:g}"
    p = lambda tag: args.outdir / f"{args.prefix}_{tag}"
    original = args.outdir / f"{stem}_hg19.bedpe"
    converted = args.outdir / f"{stem}_hg38.bedpe"
    cis = args.outdir / f"{args.prefix}_cis_score{args.min_score:g}_hg38.bedpe"
    bait_hg19, oe_hg19 = p("bait_hg19.bed"), p("otherEnd_hg19.bed")
    bait_hg38, oe_hg38 = p("bait_hg38.bed"), p("otherEnd_hg38.bed")
    bait_unmapped, oe_unmapped = p("bait_unmapped.bed"), p("otherEnd_unmapped.bed")
    failed = p("liftover_failed.tsv")

    selected = extract(args.matrix, original, bait_hg19, oe_hg19, args)
    run_liftover(args.liftover, bait_hg19, args.chain, bait_hg38, bait_unmapped)
    run_liftover(args.liftover, oe_hg19, args.chain, oe_hg38, oe_unmapped)
    counts, bait_coordinates = reconstruct(
        original, load_mapped(bait_hg38), load_mapped(oe_hg38), converted, cis, failed
    )
    if args.expected_hg38_bait and bait_coordinates != [args.expected_hg38_bait]:
        raise ValueError(f"Converted bait(s) {bait_coordinates} differ from expected {args.expected_hg38_bait}")

    summary = {
        "source": "Javierre et al. 2016 human promoter capture Hi-C peak matrix",
        "input_assembly": "hg19/GRCh37",
        "output_assembly": "hg38/GRCh38",
        "input_coordinate_convention": "1-based inclusive HindIII restriction fragments",
        "output_coordinate_convention": "BEDPE 0-based half-open intervals",
        "score_definition": "Original cell-type-specific CHiCAGO score, unchanged by liftover",
        "matrix": str(args.matrix),
        "chain": str(args.chain),
        "bait_id": args.bait_id,
        "bait_chrom": chrom(args.bait_chrom),
        "cell_column": args.cell_column,
        "minimum_score": args.min_score,
        "selected_interactions": selected,
        "converted_interactions": counts["converted"],
        "failed_interactions": counts["failed"],
        "converted_cis_interactions": counts["cis"],
        "hg38_bait_bed_intervals": bait_coordinates,
        "outputs": {"hg19_bedpe": str(original), "hg38_bedpe": str(converted),
                    "hg38_cis_bedpe": str(cis), "failed_report": str(failed),
                    "bait_unmapped": str(bait_unmapped), "other_end_unmapped": str(oe_unmapped)},
    }
    summary_path = p("liftover_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Selected in hg19: {selected}")
    print(f"Converted to hg38: {counts['converted']}")
    print(f"Failed: {counts['failed']}")
    print(f"Converted cis interactions: {counts['cis']}")
    print("Converted bait (BED coordinates):", ", ".join(bait_coordinates))
    print("BEDPE:", converted)
    print("Cis BEDPE:", cis)
    print("Report:", summary_path)


if __name__ == "__main__":
    main()
