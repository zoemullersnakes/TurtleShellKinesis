#TRIMAL=/home/user/Analysis/trimAl_Linux_x86-64/trimal
TREE=/media/user/Draco/turtles/protein_files/species_convergence3.nwk
#TRIMDIR=/media/user/Draco/turtles/protein_files/trimmed
OUTDIR=/media/user/Draco/turtles/protein_files/absrel
mkdir -p "$OUTDIR"
export TREE OUTDIR


run_ph() {
  in="$1"
  base=$(basename "$in" _trim.phy)

  out="$OUTDIR/${base}.absrel.json"
  [ -s "$out" ] && return 0                     # skip finished loci

  hyphy CPU=1 absrel --alignment "$in" --tree "$TREE" \
        --output "$out"  > "$OUTDIR/${base}.log" 2>&1
}
export -f run_ph

cd /media/user/Draco/turtles/protein_files/trimmed
parallel -j 40 --bar --joblog "$OUTDIR/hyphy_absrel.log" run_ph ::: *.phy

