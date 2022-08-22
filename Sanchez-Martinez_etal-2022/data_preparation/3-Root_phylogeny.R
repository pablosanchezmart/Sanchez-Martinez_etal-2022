#### ROOT PHYLOGENY ---------------------------------------------------------------------------------------------------- ####

## P Sanchez-Martinez

#### PACKAGES ------------------------------------------------------------------------------------------------------------- ####

library(ape)
library(V.PhyloMaker)

#### PHYLOGENY ----------------------------------------------------------------------------------------------------------- ####

tree1 <- read.tree("data/phylogeny/tree_final_46708sp_phylogeny.tre")

is.ultrametric(tree1)
is.rooted(tree1)
is.binary.tree(tree1)

#not necessarily expected to be binary, but should be rooted if ultrametric...

Nnode(tree1); Ntip(tree1)
#very big tree, albeit with LOTS of polytomies.

sort(tree1$node.label,dec=T)
#some weird node labels... let's assign numeric ones that should match what is in edge object

tree1$node.label[] <- seq(Ntip(tree1)+1 , Ntip(tree1)+Nnode(tree1))

#write that to file
write.tree(tree1,"data/phylogeny/Big_phylo.tre")

#opening in Figtree software (https://github.com/rambaut/figtree/releases)
#selecting to view node labels (then select labels) to understand nodes

#so... root node is 46709 as would be expected from looking at edge attribute of tree
tree1$edge[1:10,]

#we also see the phylogeny has some major issues, I think to do with how unsampled taxa were inserted as polytomies

#for example, we see that there is a Cupressaceae stuck on the base for some reason... let's remove that

tree2 <- drop.tip(tree1,"Fokienia_hodginsii")

#make sure I spelled that right
Ntip(tree1); Ntip(tree2)

#having gotten rid of the polytomy at the very base of the tree, it should now be rooted
is.rooted(tree2)
#indeed, now rooted

write.tree(tree2, "data/phylogeny/bigphylo2.tre")