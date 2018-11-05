library(ape)
library(dendextend)


tree_shared_int = read.tree('tree1_shared_introgression')
tree_nonint = read.tree('outtree_upgma_nonint_ungapped')

#tree_all = root(tree_all, 'CBS432', resolve.root=T)
#tree_nonint = root(tree_nonint, 'CBS432', resolve.root=T)

ctree_all = compute.brlen(tree_shared_int, method="Grafen")
#ctree_all = chronos(tree_shared_int, model='relaxed')
ctree_nonint = chronos(tree_nonint, model='relaxed')

entanglement(ctree_all, ctree_nonint)

dendan <- dendlist(ctree_all %>% as.dendrogram, ctree_nonint %>% as.dendrogram)

dendan = untangle(dendan, method='step2')

pdf("tanglegram_sharedint1_nonint.pdf", width=7, height=7)
tanglegram(dendan, common_subtrees_color_branches = TRUE, edge.lwd=1, lwd=2, lab.cex=.6, rank_branches=T)
dev.off()



tree_shared_int = read.tree('tree2_shared_introgression')
tree_nonint = read.tree('outtree_upgma_nonint_ungapped')

#tree_all = root(tree_all, 'CBS432', resolve.root=T)
#tree_nonint = root(tree_nonint, 'CBS432', resolve.root=T)

ctree_all = compute.brlen(tree_shared_int, method="Grafen")
#ctree_all = chronos(tree_shared_int, model='relaxed')
ctree_nonint = chronos(tree_nonint, model='relaxed')

entanglement(ctree_all, ctree_nonint)

dendan <- dendlist(ctree_all %>% as.dendrogram, ctree_nonint %>% as.dendrogram)

dendan = untangle(dendan, method='step2')

pdf("tanglegram_sharedint2_nonint.pdf", width=7, height=7)
tanglegram(dendan, common_subtrees_color_branches = TRUE, edge.lwd=1, lwd=2, lab.cex=.6, rank_branches=T)
dev.off()



tree_shared_int = read.tree('tree3_shared_introgression')
tree_nonint = read.tree('outtree_upgma_nonint_ungapped')

#tree_all = root(tree_all, 'CBS432', resolve.root=T)
#tree_nonint = root(tree_nonint, 'CBS432', resolve.root=T)

ctree_all = compute.brlen(tree_shared_int, method="Grafen")
#ctree_all = chronos(tree_shared_int, model='relaxed')
ctree_nonint = chronos(tree_nonint, model='relaxed')

entanglement(ctree_all, ctree_nonint)

dendan <- dendlist(ctree_all %>% as.dendrogram, ctree_nonint %>% as.dendrogram)

dendan = untangle(dendan, method='step2')

pdf("tanglegram_sharedint3_nonint.pdf", width=7, height=7)
tanglegram(dendan, common_subtrees_color_branches = TRUE, edge.lwd=1, lwd=2, lab.cex=.6, rank_branches=T)
dev.off()

sagsdag


tree_all = read.tree('outtree_consense_all_ungapped_bootstrap')
tree_nonint = read.tree('outtree_upgma_nonint_ungapped')

#tree_all = root(tree_all, 'CBS432', resolve.root=T)
#tree_nonint = root(tree_nonint, 'CBS432', resolve.root=T)

ctree_all = chronos(tree_all, model='relaxed')
ctree_nonint = chronos(tree_nonint, model='relaxed')

entanglement(ctree_all, ctree_nonint)

dendan <- dendlist(ctree_all %>% as.dendrogram, ctree_nonint %>% as.dendrogram)

dendan = untangle(dendan, method='step2')

pdf("tanglegram_allcon_nonint.pdf", width=7, height=7)
tanglegram(dendan, common_subtrees_color_branches = TRUE, edge.lwd=1, lwd=2, lab.cex=.6, rank_branches=T)
dev.off()

adsg


tree_all = read.tree('outtree_upgma_all_ungapped')
tree_nonint = read.tree('outtree_upgma_nonint_ungapped')

#tree_all = root(tree_all, 'CBS432', resolve.root=T)
#tree_nonint = root(tree_nonint, 'CBS432', resolve.root=T)

ctree_all = chronos(tree_all, model='relaxed')
ctree_nonint = chronos(tree_nonint, model='relaxed')

entanglement(ctree_all, ctree_nonint)

dendan <- dendlist(ctree_all %>% as.dendrogram, ctree_nonint %>% as.dendrogram)

dendan = untangle(dendan, method='step2')

pdf("tanglegram_all_nonint.pdf", width=7, height=7)
tanglegram(dendan, common_subtrees_color_branches = TRUE, edge.lwd=1, lwd=2, lab.cex=.6, hang.dendrogram=T)
dev.off()



tree_all = read.tree('outtree_upgma_all_ungapped')
tree_int = read.tree('outtree_upgma_int_ungapped')

#tree_all = root(tree_all, 'CBS432', resolve.root=T)
#tree_int = root(tree_int, 'CBS432', resolve.root=T)

ctree_all = chronos(tree_all, model='relaxed')
ctree_int = chronos(tree_int, model='relaxed')

entanglement(ctree_all, ctree_int)

dendan <- dendlist(ctree_all %>% as.dendrogram, ctree_int %>% as.dendrogram)

dendan = untangle(dendan, method='step2')

pdf("tanglegram_all_int.pdf", width=7, height=7)
tanglegram(dendan, common_subtrees_color_branches = TRUE, edge.lwd=1, lwd=2, lab.cex=.6, rank_branches=T)
dev.off()



asdgasdg



trees = read.tree('outtree_dollop')
pdf("dollop_consensus.pdf", width=7, height=7)
plot(consensus(trees), cex = .4, direction="leftwards", layout=2)
dev.off()

asdf

trees = read.tree('outtree_dollop')
tree_nonint = read.tree('outtree_upgma_nonint')

t = consensus(trees)
ctree = compute.brlen(t, method="Grafen")
ctree_nonint = chronos(tree_nonint, model='relaxed')

tanglegram(ctree, ctree_nonint)


trees = read.tree('outtree_dollop')
tree_nonint = read.tree('outtree_upgma_nonint')

t = consensus(trees, p=.5)
ctree = compute.brlen(t, method="Grafen")
ctree_nonint = chronos(tree_nonint, model='relaxed')

entanglement(ctree, ctree_nonint)

dendan <- dendlist(ctree %>% as.dendrogram, ctree_nonint %>% as.dendrogram)

dendan = untangle(dendan, method='step2')

pdf("tanglegram_dollop_nonint.pdf", width=7, height=7)
tanglegram(dendan, common_subtrees_color_branches = TRUE, edge.lwd=1, lwd=2, lab.cex=.6, rank_branches=T)
dev.off()


asdg
