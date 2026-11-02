#TRIMAL=/home/user/Analysis/trimAl_Linux_x86-64/trimal
TREE=/media/user/Draco/turtles/protein_files/species_convergence2.nwk
#TRIMDIR=/media/user/Draco/turtles/protein_files/trimmed
OUTDIR=/media/user/Draco/turtles/protein_files/bustedph
mkdir -p "$OUTDIR"
export TREE OUTDIR


run_ph() {
  in="$1"
  base=$(basename "$in" _trim.phy)

  out="$OUTDIR/${base}.bustedph.json"
  [ -s "$out" ] && return 0                     # skip finished loci

  hyphy CPU=1 busted-ph --alignment "$in" --tree "$TREE" \
        --output "$out" --branches Test > "$OUTDIR/${base}.log" 2>&1
}
export -f run_ph

cd /media/user/Draco/turtles/protein_files/trimmed
parallel -j 40 --bar --joblog "$OUTDIR/hyphy_bustedph.log" run_ph ::: *.phy

