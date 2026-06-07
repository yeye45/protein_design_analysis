#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a Foldseek DB from PDB/mmCIF files and export AA + 3Di sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--pdb-dir", required=True, type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--foldseek-bin", default="foldseek")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--keep-work-dir", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def run_command(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    printable = " ".join(cmd)
    print(f"[extract_3di] running: {printable}", file=sys.stderr)
    with log_path.open("at", encoding="utf-8") as log:
        log.write(f"\n$ {printable}\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {proc.returncode}: {printable}; see {log_path}")


def try_command(cmd: list[str], log_path: Path) -> bool:
    try:
        run_command(cmd, log_path)
        return True
    except RuntimeError as exc:
        print(f"[extract_3di] warning: {exc}", file=sys.stderr)
        return False


def normalize_id(header: str) -> str:
    return header.split()[0]


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    seq_parts: list[str] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield normalize_id(header), "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(f"FASTA sequence before header in {path} at line {line_no}")
                seq_parts.append(line)
        if header is not None:
            yield normalize_id(header), "".join(seq_parts)


def fasta_to_dict(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    for seq_id, seq in iter_fasta(path):
        if seq_id in data:
            raise ValueError(f"Duplicate FASTA ID in {path}: {seq_id}")
        data[seq_id] = seq
    return data


def count_structure_files(pdb_dir: Path) -> int:
    suffixes = {".pdb", ".ent", ".cif", ".mmcif", ".bcif"}
    return sum(1 for p in pdb_dir.rglob("*") if p.is_file() and p.suffix.lower() in suffixes)


def prepare_paths(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    if not args.pdb_dir.is_dir():
        raise FileNotFoundError(f"Structure directory does not exist: {args.pdb_dir}")
    n_structures = count_structure_files(args.pdb_dir)
    if n_structures == 0:
        raise FileNotFoundError(f"No PDB/mmCIF-like structure files found under: {args.pdb_dir}")
    print(f"[extract_3di] structure_files={n_structures}", file=sys.stderr)

    if args.out_csv.exists() and not args.overwrite:
        raise FileExistsError(f"Output CSV exists; use --overwrite: {args.out_csv}")
    if args.work_dir.exists():
        if args.overwrite:
            shutil.rmtree(args.work_dir)
        else:
            raise FileExistsError(f"Work directory exists; use --overwrite or a new --work-dir: {args.work_dir}")

    args.work_dir.mkdir(parents=True)
    return args.work_dir / "structures", args.work_dir / "amino_acid.fasta", args.work_dir / "three_di.fasta", args.work_dir / "foldseek.log"


def ensure_foldseek_available(foldseek_bin: str) -> None:
    if shutil.which(foldseek_bin) is None:
        raise FileNotFoundError(f"Cannot find Foldseek executable: {foldseek_bin}")


def main() -> int:
    args = parse_args()
    ensure_foldseek_available(args.foldseek_bin)
    db, aa_fasta, ss_fasta, log_path = prepare_paths(args)

    run_command([args.foldseek_bin, "createdb", str(args.pdb_dir), str(db), "--threads", str(args.threads)], log_path)
    run_command([args.foldseek_bin, "convert2fasta", str(db), str(aa_fasta)], log_path)

    ss_db = Path(str(db) + "_ss")
    if not ss_db.exists():
        raise FileNotFoundError(f"Foldseek 3Di database was not created: {ss_db}")

    if not try_command([args.foldseek_bin, "convert2fasta", str(ss_db), str(ss_fasta)], log_path):
        main_header = Path(str(db) + "_h")
        ss_header = Path(str(db) + "_ss_h")
        if not main_header.exists():
            raise FileNotFoundError(f"Cannot locate Foldseek header DB: {main_header}")
        run_command([args.foldseek_bin, "lndb", str(main_header), str(ss_header)], log_path)
        run_command([args.foldseek_bin, "convert2fasta", str(ss_db), str(ss_fasta)], log_path)

    aa_by_id = fasta_to_dict(aa_fasta)
    written = 0
    missing_aa = 0
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("wt", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow(["Sequence_ID", "Amino_Acid_Seq", "3Di_Seq"])
        for seq_id, three_di in iter_fasta(ss_fasta):
            aa = aa_by_id.get(seq_id)
            if aa is None:
                missing_aa += 1
                aa = ""
            writer.writerow([seq_id, aa, three_di])
            written += 1
            if written % 10000 == 0:
                print(f"[extract_3di] rows_written={written}", file=sys.stderr)

    if missing_aa:
        print(f"[extract_3di] warning: 3Di records without AA sequence: {missing_aa}", file=sys.stderr)
    print(f"[extract_3di] done rows_written={written} output={args.out_csv}", file=sys.stderr)

    if not args.keep_work_dir:
        shutil.rmtree(args.work_dir)
        print(f"[extract_3di] removed work_dir={args.work_dir}", file=sys.stderr)
    else:
        print(f"[extract_3di] kept work_dir={args.work_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
