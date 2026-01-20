
### defaults
IN_ABR="card.tsv"
OUT="wrapper_out"
SQUEEZE_DIR="SqueezeMeta/sample/results"

### parse command-line options
while getopts ":i:o:s:h" opt; do
    case "${opt}" in
        i) IN_ABR="${OPTARG}" ;;
        o) OUT="${OPTARG}" ;;
        s) SQUEEZE_DIR="${OPTARG}" ;;
        h) usage ;;
        \?)
            echo "ERROR: Invalid option -${OPTARG}" >&2
            usage
            ;;
        :)
            echo "ERROR: Option -${OPTARG} requires an argument." >&2
            usage
            ;;
    esac
done


echo " "
echo " "
echo " Usage:"
echo " -i) path to Abricate results (e.g., card.tsv)"
echo " -o) path to output"
echo " -s) path to SqueezeMeta working directory"
echo " "
echo "-----------------------------------------------------"
echo " "
echo " SqueezeMeta/Abricate current wrapper configuration:"
echo " "
echo " - i: ${IN_ABR}"
echo " - o: ${OUT}"
echo " - s: ${SQUEEZE_DIR}"
echo " "
echo "-----------------------------------------------------"
echo " "
echo " "

### check Abricate results file
if [[ ! -r "${IN_ABR}" ]]; then
    echo "ERROR: ${IN_ABR} not found. Is it .tsv? Use -h flag for help." >&2
    exit 1
fi


### check SqueezeMeta directory
#if [[ ! "${SQUEEZE_DIR}" =~ /results/?$ ]]; then
#    SQUEEZE_DIR="${SQUEEZE_DIR%/}/results"
#fi

SQUEEZE_DIR="$(printf '%s' "$SQUEEZE_DIR" | sed 's#//*#/#g')"

if [[ ! -d "${SQUEEZE_DIR}" ]]; then
    echo "ERROR: Input directory '${SQUEEZE_DIR}' not found." >&2
    exit 1
fi


### rewrite Abricate TSV by adding TAXON from SqueezeMeta fun3 wranks
OUT_FILE="${OUT}"

### ensure parent dir exists if OUT includes a path
mkdir -p "$(dirname "$OUT_FILE")" 2>/dev/null || true
: > "$OUT_FILE"   # truncate

### base directory. Just in case input contains multiple samples
SQUEEZE_BASE="$(dirname "$(dirname "$SQUEEZE_DIR")")"

awk -v FS='\t' -v OFS='\t' \
    -v out="$OUT_FILE" \
    -v squeeze_dir="$SQUEEZE_DIR" \
    -v squeeze_base="$SQUEEZE_BASE" '
function pick_taxfile(sample,   f1,f2,probe) {
    # If squeeze_dir already points to *this* sample/results:
    f1 = squeeze_dir "/06." sample ".fun3.tax.wranks"

    # Otherwise assume squeeze_base contains multiple sample dirs:
    #   squeeze_base/<sample>/results/06.<sample>.fun3.tax.wranks
    f2 = squeeze_base "/" sample "/results/06." sample ".fun3.tax.wranks"

    # Probe readability by trying to read one line (then close; we reopen later).
    if ((getline probe < f1) > 0) { close(f1); return f1 }
    close(f1)
    if ((getline probe < f2) > 0) { close(f2); return f2 }
    close(f2)
    return ""
}

# GENE -> TAXON
NR == 1 || $0 ~ /^#/ {
    line = $0
    gsub(/GENE/, "TAXON", line)
    print line >> out
    next
}

{
    # Abricate columns:
    # 1: bits[0] sample/path, 2: contig, 3: start, 4: end, 6: bits[5] (GENE col in many abricate outputs)
    sample = $1
    sub(/\/.*/, "", sample)  # first path segment

    contig = $2
    start  = $3 + 0
    end    = $4 + 0

    taxfile = pick_taxfile(sample)
    if (taxfile == "") {
        print "WARN: tax file not found for sample=" sample " (tried under " squeeze_dir " and " squeeze_base ")" > "/dev/stderr"
        print contig "\t<missing_taxfile>" > "/dev/stderr"
        next
    }

    done = 0
    while ((getline p < taxfile) > 0) {
        split(p, ps, FS)

        id = ps[1]
        n_parts = split(id, parts, "_")
        c = parts[1] "_" parts[2]
        if (c != contig) continue

        range_part = parts[n_parts]
        split(range_part, r, "-")
        rs = r[1] + 0
        re = r[2] + 0

        # tolerance of <=50 bp slack on either end
        if (rs <= start || (rs - start) <= 50) {
            if (end <= re || (end - re) <= 50) {

                n_tax = split(ps[2], taxbits, "_")
                tax = taxbits[n_tax]

                out_line = $1 OFS $2 OFS $3 OFS $4 OFS $5
                out_line = out_line OFS (tax ";" $6)

                for (i = 7; i <= NF; i++) out_line = out_line OFS $i

                print out_line >> out
                done = 1
                break
            }
        }
    }
    close(taxfile)

    if (!done) {
	# print unmatched resistance genes:
        print contig "\t" taxfile > "/dev/stderr"
    }
}
' "$IN_ABR"

