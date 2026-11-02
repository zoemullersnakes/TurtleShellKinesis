#library(devtools)
#pak::pak("nclark-lab/RERconverge")
#install_github("nclark-lab/RERconverge")
mastertreefn = "/Users/zlm9173/Desktop/Turtles/species_convergence.nwk"
alignmentfn = "/Users/zlm9173/Desktop/Turtles/rerconverge"  # same alignments IQ-TREE used
outputfn = "/Users/zlm9173/Desktop/Turtles/rerconverge.phangorn.tree.nwk"

estimatePhangornTreeAll(treefile = mastertreefn, alndir = alignmentfn,
                        output.file = outputfn, type = "AA", format = "fasta")


library(RERconverge)
treefile = "/Users/zlm9173/Desktop/Turtles/rerconverge.phangorn.tree.nwk"
mastertree_file = "/Users/zlm9173/Desktop/Turtles/species_convergence.nwk"
library(ape)
mastertree_obj = read.tree(mastertree_file)
trees = readTrees(treefile, useSpecies = NULL,
                  minSpecs = 15)

kineticforeground = c("Cuotri", "Emys_orbicularis", 
                      "Terrapene_triunguis", "Emydoidea_blandingii")
kinetic = foreground2Tree(kineticforeground, trees, clade = "all", transition = "unidirectional",
                          useSpecies=trees$masterTree$tip.label)
kRERw = getAllResiduals(trees, useSpecies=trees$masterTree$tip.label,
                        transform = "sqrt", use.weights = T, norm="scale")
saveRDS(kRERw, file = "/Users/zlm9173/Desktop/Turtles/rer.rds")

phenvKinetic = tree2Paths(kinetic, trees)

corKinetic = correlateWithBinaryPhenotype(kRERw, phenvKinetic, min.sp = 10,
                                          min.pos = 2, weighted = "auto",
                                          winsorizeRER = NULL, winsorizetrait = NULL,
                                          bootstrap = TRUE, bootn=1000, sort = TRUE)
head(corKinetic[order(corKinetic$P),])

corKinetic$stat = -log10(corKinetic$P) * sign(corKinetic$Rho)
head(corKinetic[order(corKinetic$stat, decreasing = TRUE),])
plotRers(kRERw, "OG0012054_aligned", phenv=phenvKinetic)
plotRers(kRERw, "OG0010495_aligned", phenv=phenvKinetic)

write.csv(corKinetic, file = "/Users/zlm9173/Desktop/Turtles/RERconverge_correlation_results.csv", row.names = TRUE)

hist(corKinetic$P, breaks=15, xlab="Kendall P-value",
     main="P-values for correlation between 200 genes and marine environment")
stats = getStat(corKinetic)
view(stats)
stats

n <- sum(!is.na(corKinetic$P))
sum(corKinetic$P < 0.05, na.rm = TRUE)   # observed
0.05 * n                                  # expected under the null
hist(corKinetic$P, breaks = 50)
