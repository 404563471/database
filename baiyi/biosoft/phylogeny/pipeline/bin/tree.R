suppressMessages(library(ggtree))
args <- commandArgs(T)
#test.tree <- read.tree(args[1], keep.multi = TRUE)
test.tree <- read.tree(args[1])

#test.tree <- fortify(test.tree)
#test.tree <- test.tree[order(test.tree$x),]
#tr <- read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
#plot(tr)

tree.leng <- test.tree$edge.length
test.tree$edge.length <- (tree.leng-min(tree.leng))/(max(tree.leng)-min(tree.leng))

tree.plot <- ggtree(test.tree) + geom_tiplab() + 
  geom_point(color='firebrick') +
  xlim(ifelse(min(test.tree[[1]][[2]])<=0, min(test.tree[[1]][[2]]), 0), 1.5)
ggsave(args[2], plot = tree.plot, device = "pdf", dpi=300)
