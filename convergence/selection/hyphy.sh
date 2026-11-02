TRIMAL=/home/user/Analysis/trimAl_Linux_x86-64/trimal
TREE=/media/user/Draco/turtles/protein_files/species_convergence2.nwk
TRIMDIR=/media/user/Draco/turtles/protein_files/trimmed
OUTDIR=/media/user/Draco/turtles/protein_files/busted
mkdir -p "$TRIMDIR" "$OUTDIR"
export TRIMAL TREE TRIMDIR OUTDIR

run_one() {
  in="$1"
  base=$(basename "$in" .phy)
  trimmed="$TRIMDIR/${base}_trim.phy"
  out="$OUTDIR/${base}.busted.json"
  [ -s "$out" ] && return 0                     # skip finished loci
  "$TRIMAL" -in "$in" -out "$trimmed" -phylip || return 1
  hyphy CPU=1 busted --alignment "$trimmed" --tree "$TREE" \
        --output "$out" --branches Test > "$OUTDIR/${base}.log" 2>&1
}
export -f run_one

cd /media/user/Draco/turtles/protein_files/phylip
parallel -j 40 --bar --joblog "$OUTDIR/hyphy_trimal.log" run_one ::: *.phy
