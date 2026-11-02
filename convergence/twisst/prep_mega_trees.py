import os 

def prep_mega_tree(iqtree_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    trees = {}
    input = [f for f in os.listdir(iqtree_folder) if f.endswith(".treefile")]
    
    for f in input:
      path = os.path.join(iqtree_folder, f)
      name = os.path.basename(path).split("_")[0]
      with open(path, "r") as i:
        tree = i.readlines()
        trees[name] = tree
    with open(os.path.join(output_folder, "twisst.tree.nwk"), "w") as o:
      for name, tree in trees.items():
        tree_clean = "".join(tree)
        o.write(f"{tree_clean}")

    with open(os.path.join(output_folder, "rerconverge.tree.nwk"), "w") as o:
      for name, tree in trees.items():
        tree_clean = "".join(tree)
        name_clean = "".join(name)
        o.write(f"{name_clean}\t{tree_clean}")
    return

prep_mega_tree("/media/user/Draco/turtles/protein_files/amino_aligned/headers", "/media/user/Draco/turtles/protein_files/") 
