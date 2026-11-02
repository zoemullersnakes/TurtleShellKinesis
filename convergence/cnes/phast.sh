#!/bin/bash
#SBATCH -J phastCons_setup
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH --mail-user=zoe.muller@uta.edu
#SBATCH --mail-type=END,FAIL
#
# ═══════════════════════════════════════════════════════════════════
# phastCons on the full turtle HAL, with a NEUTRAL model fitted to
# 4-fold degenerate (4d) sites, run in fixed-size reference windows.
#
# RUN ONCE:
#     cd /anvil/scratch/x-zlm9173
#     mkdir -p logs
#     sbatch phast_full.sh
#
# Setup submits everything below as a chain of dependent jobs:
#   A. extract (array)  hal2maf for each window + its 4d sites
#   B. model   (1 job)  pool 4d sites genome-wide -> phyloFit -> neutral.mod
#   C. score   (array)  phastCons on each window with neutral.mod
#   D. merge   (1 job)  per-chromosome wig + genome-wide conserved BED
#
#   sbatch phast_full.sh retry   reruns unfinished windows (step C) with
#                                more memory, then merges again.
# Re-submitting is safe: finished windows are skipped.
# ═══════════════════════════════════════════════════════════════════

# ─── SETTINGS ──────────────────────────────────────────────────────
SCR=/anvil/scratch/x-zlm9173
SCRIPT=$SCR/phast.sh                 # where THIS file lives
IMG=/apps/biocontainers/images/quay.io_comparative-genomics-toolkit_cactus:v2.6.5.sif
HAL=$SCR/turtles.hal
REF=Chrysemys_picta_bellii
# Reference annotation (RefSeq GTF, .gtf or .gtf.gz). Column 1 must use the
# same sequence names as the Chrysemys genome inside the HAL (e.g. NC_...).
GTF=$SCR/Chrysemys_picta_bellii.gtf

MIN_SCAFFOLD=100000         # skip reference scaffolds shorter than this (bp)
WINDOW=10000000             # window size (bp)
MAX_PARALLEL=50             # array tasks running at once
KEEP_MAF=1                  # keep per-window MAFs after scoring

PC_ARGS="--target-coverage 0.3 --expected-length 45 --rho 0.3"

WORKER_MEM=16G
WORKER_TIME=12:00:00
RETRY_PARTITION=shared
RETRY_CPUS=16
RETRY_MEM=32G
RETRY_TIME=24:00:00

# Everything below lives in its own folder, separate from the partial run
RUN=$SCR/phast_turtles
TREE=$RUN/turtles.noanc.nh
SIZES=$RUN/ref.big.sizes
WINDOWS=$RUN/windows.txt
GTFDIR=$RUN/cds_by_chrom
MAFDIR=$RUN/maf_windows
FDDIR=$RUN/4d_windows
MOD=$RUN/neutral.mod
OUTDIR=$RUN/windows_out
FINAL=$RUN/final
# ───────────────────────────────────────────────────────────────────

set -o pipefail
mkdir -p $SCR/logs $RUN $GTFDIR $MAFDIR $FDDIR $OUTDIR $FINAL
cd $SCR
hal() { singularity exec -B $SCR $IMG "$@"; }

module --force purge
ml biocontainers phast

MODE=${1:-setup}

window_line() { awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n' "$1"; }

# Comma-separated species (tips of the tree) — the only rows allowed in a MAF
species_list() { tr -d ' \n;()' < "$TREE" | sed -E 's/:[0-9.eE+-]+//g'; }

# Clean each MAF once: remove comment lines (hal2maf writes "# hal <tree>",
# which PHAST misreads as sequence names) and any row whose genome is not a
# tree tip (e.g. Cactus ancestors such as
# "Anc11", which hal2maf can emit despite --noAncestors). Runs once per MAF.
filter_maf() {   # $1 maf
    local maf=$1
    [ -f "$maf.filtered" ] && return 0
    awk -v keep="$(species_list)" '
        BEGIN { n = split(keep, a, ","); for (i = 1; i <= n; i++) k[a[i]] = 1 }
        /^#/ && !/^##maf/ { next }      # e.g. "# hal <tree>": PHAST misreads it as names
        $1 ~ /^[sieq]$/ { g = $2; sub(/\..*/, "", g); if (!(g in k)) { drop[g]++; next } }
        { print }
        END { for (g in drop) printf "  removed %d rows of non-tip genome %s\n", drop[g], g > "/dev/stderr" }
    ' "$maf" > "$maf.tmp" || { rm -f "$maf.tmp"; return 1; }
    mv "$maf.tmp" "$maf" && touch "$maf.filtered"
}

make_maf() {   # $1 id $2 chr $3 start $4 len
    local maf=$MAFDIR/$1.maf
    if [ ! -s "$maf" ]; then
        hal hal2maf "$HAL" "$maf.tmp" --refGenome "$REF" --refSequence "$2" --start "$3" --length "$4" --unique --noAncestors --noDupes || { rm -f "$maf.tmp"; return 1; }
        mv "$maf.tmp" "$maf"
    fi
    filter_maf "$maf"
}

check_array_size() {
    local n maxarr
    n=$(wc -l < "$1")
    maxarr=$(scontrol show config | awk '/^MaxArraySize/ {print $3}')
    if [ -n "$maxarr" ] && [ "$n" -ge "$maxarr" ]; then
        echo "ERROR: $n tasks exceeds MaxArraySize=$maxarr; raise WINDOW or MIN_SCAFFOLD"; exit 1
    fi
}

# ═══════════════════════ SETUP ═════════════════════════════════════
if [ "$MODE" = setup ]; then
    echo "== SETUP $(date)"
    [ -f "$SCRIPT" ] || { echo "ERROR: set SCRIPT= to this file's path ($SCRIPT not found)"; exit 1; }
    [ -f "$HAL" ]    || { echo "ERROR: $HAL not found"; exit 1; }

    if [ ! -s "$TREE" ]; then
        hal halStats --tree "$HAL" | sed -E 's/\)Anc[0-9]+/)/g' > "$TREE" || exit 1
    fi
    grep -q ';' "$TREE" || { echo "ERROR: tree in $TREE is empty or invalid"; rm -f "$TREE"; exit 1; }
    NSP=$(tr ',' '\n' < "$TREE" | wc -l)
    echo "Tree ($NSP species): $(cat $TREE)"

    hal halStats --chromSizes "$REF" "$HAL" | awk -v m=$MIN_SCAFFOLD '$2 >= m' > "$SIZES" || exit 1
    [ -s "$SIZES" ] || { echo "ERROR: no $REF scaffolds >= $MIN_SCAFFOLD bp"; exit 1; }

    awk -v w=$WINDOW '{ for (s = 0; s < $2; s += w) {
            len = ($2 - s < w) ? $2 - s : w
            printf "%s__%010d\t%s\t%d\t%d\n", $1, s, $1, s, len } }' "$SIZES" > "$WINDOWS"
    echo "Scaffolds: $(wc -l < $SIZES)   windows: $(wc -l < $WINDOWS) (${WINDOW} bp)"
    check_array_size "$WINDOWS"

    # Neutral model inputs: CDS lines per chromosome, column 1 renamed to the
    # species (PHAST matches features to sequences by species name)
    if [ ! -s "$MOD" ]; then
        [ -s "$GTF" ] || { echo "ERROR: GTF not found: $GTF"; exit 1; }
        rm -f $GTFDIR/*.gtf
        CAT=cat; case "$GTF" in *.gz) CAT=zcat;; esac
        $CAT "$GTF" | awk -F'\t' -v OFS='\t' -v sp="$REF" -v d="$GTFDIR" -v sizes="$SIZES" '
            BEGIN { while ((getline l < sizes) > 0) { split(l, a, "\t"); keep[a[1]] = 1 } }
            $3 == "CDS" && ($1 in keep) { c = $1; $1 = sp; print > (d "/" c ".gtf") }'
        NCHR=$(ls $GTFDIR/*.gtf 2>/dev/null | wc -l)
        echo "Chromosomes with CDS annotation: $NCHR / $(wc -l < $SIZES)"
        if [ "$NCHR" -eq 0 ]; then
            echo "ERROR: no GTF CDS lines match the HAL's $REF sequence names."
            echo "  HAL names: $(head -3 $SIZES | cut -f1 | tr '\n' ' ')"
            echo "  GTF names: $($CAT "$GTF" | grep -v '^#' | cut -f1 | sort -u | head -3 | tr '\n' ' ')"
            exit 1
        fi
    fi

    N=$(wc -l < "$WINDOWS")
    DEP=""
    if [ ! -s "$MOD" ]; then
        A=$(sbatch --parsable -J phast_extract --array=1-${N}%${MAX_PARALLEL} -t $WORKER_TIME --mem=$WORKER_MEM \
            -o $SCR/logs/phast_extract_%A_%a.out --mail-type=NONE "$SCRIPT" extract) || exit 1
        B=$(sbatch --parsable -J phast_model --dependency=afterany:$A -t 06:00:00 --mem=32G \
            -o $SCR/logs/phast_model_%j.out "$SCRIPT" model) || exit 1
        DEP="--dependency=afterok:$B"
        echo "Submitted extract array $A, model job $B"
    else
        echo "Using existing model $MOD"
    fi
    C=$(sbatch --parsable -J phast_score --array=1-${N}%${MAX_PARALLEL} $DEP -t $WORKER_TIME --mem=$WORKER_MEM \
        -o $SCR/logs/phast_score_%A_%a.out --mail-type=NONE "$SCRIPT" score "$WINDOWS") || exit 1
    D=$(sbatch --parsable -J phast_merge --dependency=afterany:$C -t 02:00:00 --mem=8G \
        -o $SCR/logs/phast_merge_%j.out "$SCRIPT" merge) || exit 1
    echo "Submitted score array $C and merge job $D"
    exit 0
fi

# ═══════════════════════ A. EXTRACT ════════════════════════════════
if [ "$MODE" = extract ]; then
    read -r ID CHR START LEN < <(window_line "$WINDOWS")
    [ -n "$ID" ] || { echo "ERROR: no window on line $SLURM_ARRAY_TASK_ID"; exit 1; }
    echo "== Extract $SLURM_ARRAY_TASK_ID -> $CHR:$START+$LEN  $(date)"
    make_maf "$ID" "$CHR" "$START" "$LEN" || exit 1
    echo "MAF blocks: $(grep -c '^a' $MAFDIR/$ID.maf)"

    FEAT=$GTFDIR/$CHR.gtf
    OUT=$FDDIR/$ID.4d.ss
    if [ ! -s "$FEAT" ] || [ "$(grep -c '^a' $MAFDIR/$ID.maf)" -eq 0 ]; then
        echo "No CDS or no alignment in this window; no 4d sites."; : > "$OUT"; exit 0
    fi
    # only CDS overlapping this window (1-based GTF vs 0-based window start)
    awk -v s=$START -v e=$((START + LEN)) -F'\t' '$5 > s && $4 <= e' "$FEAT" > $FDDIR/$ID.gtf
    if [ ! -s $FDDIR/$ID.gtf ]; then
        echo "No CDS in this window."; : > "$OUT"; rm -f $FDDIR/$ID.gtf; exit 0
    fi
    msa_view $MAFDIR/$ID.maf --in-format MAF --4d --features $FDDIR/$ID.gtf > $FDDIR/$ID.4dc.ss || exit 1
    msa_view $FDDIR/$ID.4dc.ss --in-format SS --out-format SS --tuple-size 1 > "$OUT" || exit 1
    rm -f $FDDIR/$ID.4dc.ss $FDDIR/$ID.gtf
    echo "4d sites: $(awk '/^LENGTH/ {print $3}' $OUT)"
    exit 0
fi

# ═══════════════════════ B. MODEL ══════════════════════════════════
if [ "$MODE" = model ]; then
    echo "== MODEL $(date)"
    SPECIES=$(species_list)
    # keep files with sites whose sequence names are all tree tips
    LIST=$(for f in $FDDIR/*.4d.ss; do
               [ -s "$f" ] || continue
               awk '/^LENGTH/ && $3 > 0 {found=1} END {exit !found}' "$f" || continue
               bad=$(awk -v keep="$SPECIES" -F' = ' '/^NAMES/ {
                         n = split(keep, a, ","); for (i = 1; i <= n; i++) k[a[i]] = 1
                         m = split($2, s, ","); for (i = 1; i <= m; i++) if (!(s[i] in k)) printf "%s ", s[i] }' "$f")
               if [ -n "$bad" ]; then echo "SKIPPING $(basename $f): non-tip names: $bad" >&2; continue; fi
               echo "$f"
           done)
    NF=$(echo "$LIST" | grep -c .)
    echo "Windows with 4d sites: $NF / $(wc -l < $WINDOWS)"
    [ "$NF" -gt 0 ] || { echo "ERROR: no 4d sites extracted; check logs/phast_extract_*"; exit 1; }
    msa_view --aggregate "$SPECIES" --in-format SS --out-format SS --unordered-ss $LIST > $RUN/all.4d.ss || exit 1
    echo "Total 4d sites: $(awk '/^LENGTH/ {print $3}' $RUN/all.4d.ss)"
    phyloFit --tree "$TREE" --msa-format SS --subst-mod REV --out-root "${MOD%.mod}" $RUN/all.4d.ss || exit 1
    echo "Neutral model:"; grep -E "^(BACKGROUND|TREE)" "$MOD"
    exit 0
fi

# ═══════════════════════ C. SCORE ══════════════════════════════════
if [ "$MODE" = score ]; then
    LIST=${2:-$WINDOWS}
    read -r ID CHR START LEN < <(window_line "$LIST")
    [ -n "$ID" ] || { echo "ERROR: no window on line $SLURM_ARRAY_TASK_ID of $LIST"; exit 1; }
    WIG=$OUTDIR/$ID.wig; BED=$OUTDIR/$ID.bed
    echo "== Score $SLURM_ARRAY_TASK_ID -> $CHR:$START+$LEN  $(date)"
    [ -f "$OUTDIR/$ID.done" ] && { echo "Already finished; skipping."; exit 0; }
    [ -s "$MOD" ] || { echo "ERROR: model $MOD missing"; exit 1; }
    make_maf "$ID" "$CHR" "$START" "$LEN" || exit 1

    if [ "$(grep -c '^a' $MAFDIR/$ID.maf)" -eq 0 ]; then
        echo "No alignment in this window."; : > "$WIG"; : > "$BED"; touch "$OUTDIR/$ID.done"; exit 0
    fi

    phastCons --msa-format MAF $PC_ARGS --seqname "$CHR" --most-conserved "$BED" "$MAFDIR/$ID.maf" "$MOD" > "$WIG" || exit 1

    # guard: coordinates must be chromosome-wide
    WSTART=$(awk '/^fixedStep/ { for (i = 1; i <= NF; i++) if ($i ~ /^start=/) { print substr($i, 7); exit } }' "$WIG")
    if [ "$START" -gt 0 ] && [ -n "$WSTART" ] && [ "$WSTART" -le "$START" ]; then
        echo "NOTE: window-relative coordinates detected; shifting by $START"
        awk -v o=$START '/^fixedStep/ { for (i = 1; i <= NF; i++) if ($i ~ /^start=/) $i = "start=" substr($i, 7) + o } { print }' "$WIG" > "$WIG.tmp" && mv "$WIG.tmp" "$WIG"
        awk -v o=$START 'BEGIN{OFS="\t"} { $2 += o; $3 += o; print }' "$BED" > "$BED.tmp" && mv "$BED.tmp" "$BED"
    fi
    echo "Conserved elements in window: $(wc -l < "$BED")"
    touch "$OUTDIR/$ID.done"
    [ "$KEEP_MAF" = 0 ] && rm -f "$MAFDIR/$ID.maf"
    echo "Done $(date)"
    exit 0
fi

# ═══════════════════════ RETRY ═════════════════════════════════════
if [ "$MODE" = retry ]; then
    echo "== RETRY $(date)"
    LIST=$RUN/retry.txt
    awk -v d="$OUTDIR" '{ f = d "/" $1 ".done"; if (system("[ -f \"" f "\" ]") != 0) print }' "$WINDOWS" > "$LIST"
    [ -s "$LIST" ] || { echo "Nothing to retry: all windows finished."; exit 0; }
    [ -s "$MOD" ] || { echo "ERROR: $MOD missing; the model step must succeed first (see logs/phast_model_*)"; exit 1; }
    N=$(wc -l < "$LIST")
    echo "Retrying $N window(s) with $RETRY_CPUS CPUs / $RETRY_MEM on $RETRY_PARTITION"
    C=$(sbatch --parsable -J phast_retry --array=1-${N}%${MAX_PARALLEL} -p $RETRY_PARTITION \
        --cpus-per-task=$RETRY_CPUS --mem=$RETRY_MEM -t $RETRY_TIME \
        -o $SCR/logs/phast_retry_%A_%a.out --mail-type=NONE "$SCRIPT" score "$LIST") || exit 1
    D=$(sbatch --parsable -J phast_merge --dependency=afterany:$C -t 02:00:00 --mem=8G \
        -o $SCR/logs/phast_merge_%j.out "$SCRIPT" merge) || exit 1
    echo "Submitted retry array $C and merge job $D"
    exit 0
fi

# ═══════════════════════ D. MERGE ══════════════════════════════════
if [ "$MODE" = merge ]; then
    echo "== MERGE $(date)"
    NW=$(wc -l < "$WINDOWS")
    awk -v d="$OUTDIR" '{ f = d "/" $1 ".done"; if (system("[ -f \"" f "\" ]") != 0) print $1 }' "$WINDOWS" > $RUN/failed.txt
    NFAIL=$(wc -l < $RUN/failed.txt)
    echo "Windows finished: $((NW - NFAIL)) / $NW"
    [ "$NFAIL" -gt 0 ] && echo "Unfinished windows in $RUN/failed.txt; rerun with: sbatch $SCRIPT retry"

    rm -f $FINAL/*.wig
    while read -r ID CHR START LEN; do
        [ -s "$OUTDIR/$ID.wig" ] && cat "$OUTDIR/$ID.wig" >> "$FINAL/$CHR.wig"
    done < "$WINDOWS"

    cat $OUTDIR/*.bed 2>/dev/null | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}
        { if ($1 == c && $2 <= e) { if ($3 > e) e = $3; if ($5 > sc) sc = $5; next }
          if (c != "") print c, s, e, c ":" s "-" e, sc
          c = $1; s = $2; e = $3; sc = $5 }
        END { if (c != "") print c, s, e, c ":" s "-" e, sc }' > $FINAL/most-cons.all.bed

    echo "Conserved elements: $(wc -l < $FINAL/most-cons.all.bed) -> $FINAL/most-cons.all.bed"
    awk -v g=$(awk '{s += $2} END {print s}' $SIZES) '{s += $3 - $2}
        END {printf "Total conserved bp: %d (%.2f%% of analysed reference)\n", s, 100 * s / g}' $FINAL/most-cons.all.bed
    echo "Per-chromosome wig files: $FINAL/*.wig"
    exit 0
fi

echo "Unknown mode: $MODE  (use: setup | retry)"; exit 1