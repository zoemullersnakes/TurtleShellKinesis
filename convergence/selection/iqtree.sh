run_iqtree() {
  in="$1"
  base=$(basename "$in" _AA.fa)
  iqtree -s $in -m MFP -B 1000
}

export -f run_iqtree

cd '/media/user/Draco/turtles/protein_files/amino_aligned/headers'
parallel -j 60 --bar --joblog iqtree.log run_iqtree ::: *_AA.fa
