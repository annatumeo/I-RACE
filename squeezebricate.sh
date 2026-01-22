#!/usr/bin/env bash
set -euo pipefail

### defaults
IN_ABR="card.tsv"
OUT="squeezebricate_out.tsv"
SQUEEZE_DIR="."
UNMATCHED_FILE="unmatched.tsv"

usage() {
  cat >&2 <<'EOF'

Welcome to SqueezeBricate :)

This wrapper adds a TAXON column (derived from SqueezeMeta) to Abricate TSV output.

Usage: squeezebricate_cached.sh -i <abricate.tsv> -o <out.tsv> -s </path/to/your/SqueezeMeta/working/dir> [-u <unmatched_file>]

Options:
  -i  Abricate results TSV (default: card.tsv)
  -o  Output TSV (default: squeezebricate_out.tsv)
  -s  Path to SqueezeMeta working directory (default: ".")
  -u  File to write unmatched hits (default: unmatched)
  -h  Show help

EOF
  exit 2
}

while getopts ":i:o:s:u:h" opt; do
  case "${opt}" in
    i) IN_ABR="${OPTARG}" ;;
    o) OUT="${OPTARG}" ;;
    s) SQUEEZE_DIR="${OPTARG}" ;;
    u) UNMATCHED_FILE="${OPTARG}" ;;
    h) usage ;;
    \?) echo "ERROR: Invalid option -${OPTARG}" >&2; usage ;;
    :)  echo "ERROR: Option -${OPTARG} requires an argument." >&2; usage ;;
  esac
done

if [[ ! -r "${IN_ABR}" ]]; then
  echo "ERROR: ${IN_ABR} not found or not readable." >&2
  exit 1
fi

SQUEEZE_DIR="$(printf '%s' "$SQUEEZE_DIR" | sed 's#//*#/#g')"
if [[ ! -d "${SQUEEZE_DIR}" ]]; then
  echo "ERROR: Input directory '${SQUEEZE_DIR}' not found." >&2
  exit 1
fi

command -v gawk >/dev/null 2>&1 || {
  echo "ERROR: gawk is required for this optimized wrapper." >&2
  exit 1
}

mkdir -p "$(dirname "$OUT")" 2>/dev/null || true
: > "$OUT"
: > "$UNMATCHED_FILE"

SQUEEZE_BASE="$(dirname "$(dirname "$SQUEEZE_DIR")")"

# GENE -> TAXON
awk -v FS='\t' -v OFS='\t' '
  NR==1 || $0 ~ /^#/ {
    line=$0
    gsub(/GENE/,"TAXON",line)
    print line
    next
  }
  { exit }
' "$IN_ABR" >> "$OUT"

printf "SAMPLE\tCONTIG\tSTART\tEND\tGENE\tREASON\tORIGINAL_ROW\n" >> "$UNMATCHED_FILE"

TMP_SORTED="$(mktemp)"
awk -v FS='\t' -v OFS='\t' '
  NR==1 || $0 ~ /^#/ { next }
  {
    s=$1
    sub(/\/.*/,"",s)
    print s, $0
  }
' "$IN_ABR" | sort -t $'\t' -k1,1 > "$TMP_SORTED"

gawk -v FS='\t' -v OFS='\t' \
     -v out="$OUT" \
     -v unmatched="$UNMATCHED_FILE" \
     -v squeeze_dir="$SQUEEZE_DIR" \
     -v squeeze_base="$SQUEEZE_BASE" '

function pick_taxfile(sample,f1,f2,probe) {
    f1 = squeeze_dir "/06." sample ".fun3.tax.wranks"
    f2 = squeeze_base "/" sample "/results/06." sample ".fun3.tax.wranks"

    if ((getline probe < f1) > 0) { close(f1); return f1 }
    close(f1)
    if ((getline probe < f2) > 0) { close(f2); return f2 }
    close(f2)
    return ""
}

function reset_cache() {
    delete cnt
    delete st
    delete en
    delete tx
}

function load_taxfile(sample,taxfile,p,ps,id,parts,n_parts,c,range_part,r,rs,re,n_tax,taxbits,tax) {
    taxfile = pick_taxfile(sample)
    if (taxfile == "") {
        # Warn once per sample
        print "WARN: tax file not found for sample=" sample " (tried under " squeeze_dir " and " squeeze_base ")" > "/dev/stderr"
        return 0
    }

    while ((getline p < taxfile) > 0) {
        split(p, ps, FS)
        id = ps[1]

        # contig parsing matches original: first two '_' parts
        n_parts = split(id, parts, "_")
        if (n_parts < 3) continue
        c = parts[1] "_" parts[2]

        # range is last '_' part like "123-456"
        range_part = parts[n_parts]
        split(range_part, r, "-")
        rs = r[1] + 0
        re = r[2] + 0

        # taxon is last '_' part of column 2
        n_tax = split(ps[2], taxbits, "_")
        tax = taxbits[n_tax]

        cnt[c]++
        st[c, cnt[c]] = rs
        en[c, cnt[c]] = re
        tx[c, cnt[c]] = tax
    }
    close(taxfile)
    return 1
}

function write_unmatched(sample, contig, start, end, gene, reason,   orig, j) {
    orig = $2
    for (j = 3; j <= NF; j++) orig = orig OFS $j
    print sample, contig, start, end, gene, reason, orig >> unmatched
}

BEGIN {
    cur_sample = ""
    loaded = 0
    n_missing_tax = 0
    n_no_match = 0
}

{

    sample_key = $1

    if (sample_key != cur_sample) {
        reset_cache()
        loaded = load_taxfile(sample_key)
        cur_sample = sample_key
    }

    contig = $3
    start  = $4 + 0
    end    = $5 + 0
    gene   = $7   # original $6

    if (!loaded) {
        n_missing_tax++
        write_unmatched(cur_sample, contig, start, end, gene, "missing_taxfile")
        next
    }

    # tolerance is <=50 bp on either end
    done = 0
    n = cnt[contig]
    for (i = 1; i <= n; i++) {
        rs = st[contig, i]
        re = en[contig, i]
        if (rs <= start || (rs - start) <= 50) {
            if (end <= re || (end - re) <= 50) {
                tax = tx[contig, i]

                out_line = $2 OFS $3 OFS $4 OFS $5 OFS $6
                out_line = out_line OFS (tax ";" $7)
                for (f = 8; f <= NF; f++) out_line = out_line OFS $f

                print out_line >> out
                done = 1
                break
            }
        }
    }

    if (!done) {
        n_no_match++
        write_unmatched(cur_sample, contig, start, end, gene, "no_match")
    }
}

END {
    print "INFO: unmatched (missing_taxfile)=" n_missing_tax ", unmatched (no_match)=" n_no_match > "/dev/stderr"
}
' "$TMP_SORTED"

rm -f "$TMP_SORTED"

echo "Wrote: $OUT" >&2
echo "Wrote: $UNMATCHED_FILE" >&2

