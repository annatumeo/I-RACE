#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import subprocess
import sys

# ---- SETTINGS ----
CPUS = 24

MIN_COMPLETENESS = 94.0
MAX_CONTAMINATION = 5.0
# ---------------------


def safe_float(x: str) -> float | None:
    x = (x or "").strip()
    if not x:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def has_species_level(consensus_field: str) -> bool:
    return consensus_field.strip().startswith("s_")


def extract_consensus_field(tax_path: Path) -> str | None:
    try:
        lines = tax_path.read_text().splitlines()
    except FileNotFoundError:
        print(f"[WARN] Missing tax file: {tax_path}", file=sys.stderr)
        return None

    if not lines:
        print(f"[WARN] Empty tax file: {tax_path}", file=sys.stderr)
        return None

    parts = lines[-1].split(";")
    if len(parts) < 2:
        print(f"[WARN] Unexpected tax format (no ';' fields): {tax_path}", file=sys.stderr)
        return None

    return parts[-2].strip()


def run_cmd(cmd: list[str]) -> None:
    print("[CMD]", " ".join(cmd))
    subprocess.run(cmd, check=True)


def read_bintable_header_and_rows(bintable_path: Path):
    with bintable_path.open("r", encoding="utf-8", errors="replace") as fh:
        # 1) Skip first commented line (very first line)
        first = fh.readline()
        if not first:
            return None, []
        # If it *isn't* commented for some reason, we still continue; not fatal.

        # 2) Read header line
        header_line = fh.readline()
        if not header_line:
            return None, []

        header = [h.strip() for h in header_line.rstrip("\n").split("\t")]
        rows = []
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            # allow occasional comment lines elsewhere too
            if line.startswith("#"):
                continue
            rows.append(line.split("\t"))
        return header, rows


def main() -> None:
    base = Path.cwd()

    projects_raw = input("Enter a list of projects (comma-separated): ").strip()
    projects = [p.strip() for p in projects_raw.split(",") if p.strip()]
    if not projects:
        print("No projects provided. Exiting.", file=sys.stderr)
        sys.exit(1)

    for p in projects:
        results_dir = base / p / "results"
        bins_dir = results_dir / "bins"
        bintable_path = results_dir / f"18.{p}.bintable"

        if not bintable_path.exists():
            print(f"[WARN] Bintable not found for project '{p}': {bintable_path}", file=sys.stderr)
            continue

        # Ensure directory structure exists
        hq_dir = bins_dir / "hq"
        unclassified_dir = hq_dir / "unclassified"
        genomes_dir = unclassified_dir / "gtdbtk_genomes_dir"
        out_dir = unclassified_dir / "gtdbtk_results"
        for d in (hq_dir, unclassified_dir, genomes_dir, out_dir):
            d.mkdir(parents=True, exist_ok=True)

        header, rows = read_bintable_header_and_rows(bintable_path)
        if not header:
            print(f"[WARN] Could not read header/rows from: {bintable_path}", file=sys.stderr)
            continue

        # Find required columns by name
        # (strip + exact match; if your header varies in case, we can lowercase-match)
        try:
            idx_bin = header.index(header[0])  # first column is bin id in your original
            idx_comp = header.index("Completeness")
            idx_cont = header.index("Contamination")
        except ValueError as e:
            print(f"[ERROR] Missing expected column in header for {bintable_path}: {e}", file=sys.stderr)
            print(f"[INFO] Header columns: {header}", file=sys.stderr)
            continue

        genomes_added = 0

        for cols in rows:
            # Basic row sanity
            if len(cols) <= max(idx_comp, idx_cont, idx_bin):
                print(f"[WARN] Skipping short row in {bintable_path}: {cols}", file=sys.stderr)
                continue

            row_text = "\t".join(cols)
            if ("concoct" not in row_text) and ("metabat2" not in row_text):
                continue

            bin_id = cols[idx_bin].strip()
            comp = safe_float(cols[idx_comp])
            cont = safe_float(cols[idx_cont])
            if comp is None or cont is None:
                continue

            if comp <= MIN_COMPLETENESS or cont >= MAX_CONTAMINATION:
                continue

            tax_file = bins_dir / f"{bin_id}.fa.tax"
            fasta_file = bins_dir / f"{bin_id}.fa"

            # Copy HQ tax file
            if tax_file.exists():
                run_cmd(["cp", str(tax_file), str(hq_dir)])
            else:
                print(f"[WARN] HQ candidate missing tax file: {tax_file}", file=sys.stderr)

            consensus = extract_consensus_field(tax_file)
            species_assigned = has_species_level(consensus) if consensus is not None else False

            if not species_assigned:
                if tax_file.exists():
                    run_cmd(["cp", str(tax_file), str(unclassified_dir)])
                if fasta_file.exists():
                    dest = genomes_dir / f"{bin_id}.fna"
                    run_cmd(["cp", str(fasta_file), str(dest)])
                    genomes_added += 1
                else:
                    print(f"[WARN] Missing fasta for unclassified bin: {fasta_file}", file=sys.stderr)

        if genomes_added == 0:
            print(f"[INFO] Project '{p}': no unclassified HQ genomes found; skipping GTDB-Tk.")
            continue

        cmd = [
        "gtdbtk", "classify_wf",
        "--genome_dir", str(genomes_dir),
        "--out_dir", str(out_dir),
        "--cpus", str(CPUS),
	]

        run_cmd(cmd)

    print("[DONE]")


if __name__ == "__main__":
    main()

